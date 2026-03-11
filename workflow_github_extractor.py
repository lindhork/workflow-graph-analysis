#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import requests
from github import Github
from github.GithubException import GithubException


ENGINE = "nextflow"  # "snakemake", "nextflow", or "both"
MAX_REPOS = 100
MIN_STARS = 5
MONTHS_ACTIVE = 24
COMMITS_PER_REPO = 100
ISSUES_PER_REPO = 50
PRS_PER_REPO = 50
MAX_WORKERS = 6
CLONE_DEPTH = 100
REQUEST_TIMEOUT = 60
SLEEP_BETWEEN_GRAPHQL_CALLS = 0.2


UPDATE_EXISTING_REPOS = False

ENGINE = ENGINE.lower()
if ENGINE not in {"snakemake", "nextflow", "both"}:
    raise ValueError("ENGINE must be 'snakemake', 'nextflow' or 'both'")

BASE_REPOS_DIR = Path("cloned_repos")
BASE_DATA_DIR = Path("data")
BASE_RAW_DIR = Path("raw")
LOG_DIR = Path("logs")

GITHUB_GRAPHQL = "https://api.github.com/graphql"
TOKEN = os.getenv("GITHUB_TOKEN")

if not TOKEN:
    raise RuntimeError("Please set GITHUB_TOKEN as an environment variable.")

gh = Github(TOKEN, per_page=100)


RULE_RX = re.compile(r"(?ms)^rule\s+([^\s:]+)\s*:(.*?)(?=^rule\s|^checkpoint\s|\Z)")
CHECKPOINT_RX = re.compile(r"(?ms)^checkpoint\s+([^\s:]+)\s*:(.*?)(?=^rule\s|^checkpoint\s|\Z)")
PROCESS_RX = re.compile(r"(?ms)^process\s+([^\s\{]+)\s*\{(.*?)(?=^process\s|^workflow\s|\Z)")
WORKFLOW_BLOCK_RX = re.compile(r"(?ms)^workflow\s*(?:\w+)?\s*\{(.*?)\}")

SECTION_INPUT_RX = re.compile(r"(?ms)^\s*input\s*:\s*(.*?)(?=^\s*[a-zA-Z_]+\s*:|\Z)")
SECTION_OUTPUT_RX = re.compile(r"(?ms)^\s*output\s*:\s*(.*?)(?=^\s*[a-zA-Z_]+\s*:|\Z)")
RULE_REF_RX = re.compile(r"rules\.([A-Za-z_][A-Za-z0-9_]*)\.output")
QUOTED_STRING_RX = re.compile(r"['\"]([^'\"]+)['\"]")
CHANNEL_NAME_RX = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\b")

SNAKEMAKE_REPO_PATTERNS = ("Snakefile", ".smk")
SNAKEMAKE_TASK_PATTERNS = ("Snakefile", ".smk")

NEXTFLOW_REPO_PATTERNS = ("main.nf", ".nf", "nextflow.config")
NEXTFLOW_TASK_PATTERNS = ("main.nf", ".nf")


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def get_engine_dirs(engine: str) -> Tuple[Path, Path, Path]:
    return (
        BASE_REPOS_DIR / engine,
        BASE_DATA_DIR / engine,
        BASE_RAW_DIR / engine,
    )


def ensure_dirs(engine: str) -> Tuple[Path, Path, Path]:
    repos_dir, data_dir, raw_dir = get_engine_dirs(engine)
    for d in (repos_dir, data_dir, raw_dir, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)
    return repos_dir, data_dir, raw_dir


def sha1_text(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8", errors="ignore")).hexdigest()


def stable_id(*parts: object) -> str:
    joined = "::".join(str(p) for p in parts)
    return sha1_text(joined)


def safe_repo_dir_name(full_name: str) -> str:
    return full_name.replace("/", "__")


def read_text_file(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def log_message(message: str) -> None:
    print(message)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with (LOG_DIR / "extractor.log").open("a", encoding="utf-8") as fh:
        fh.write(f"{datetime.now().isoformat()} {message}\n")


def safe_get(d: object, *keys: str, default=None):
    cur = d
    for key in keys:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(key)
        if cur is None:
            return default
    return cur


def safe_nodes(d: object, *keys: str) -> List[dict]:
    value = safe_get(d, *keys, default=[])
    return value if isinstance(value, list) else []


def sanitize_for_csv(text: str) -> str:
    text = " ".join(text.strip().split())
    text = text.replace('"', "'")
    text = text.replace("\r", " ")
    text = text.replace("\x00", "")
    return text


def matches_patterns(path: str, patterns: Sequence[str]) -> bool:
    name = Path(path).name
    lower_path = path.lower()
    lower_name = name.lower()

    for pattern in patterns:
        pattern_lower = pattern.lower()
        if pattern.startswith("."):
            if lower_path.endswith(pattern_lower):
                return True
        else:
            if lower_name == pattern_lower:
                return True
    return False


def matches_workflow_file(path: str, engine: str, purpose: str = "repo") -> bool:
    if engine == "snakemake":
        patterns = SNAKEMAKE_REPO_PATTERNS if purpose == "repo" else SNAKEMAKE_TASK_PATTERNS
        return matches_patterns(path, patterns)

    if engine == "nextflow":
        patterns = NEXTFLOW_REPO_PATTERNS if purpose == "repo" else NEXTFLOW_TASK_PATTERNS
        return matches_patterns(path, patterns)

    return False


@dataclass
class RepoCandidate:
    repo_id: str
    full_name: str
    name: str
    owner: str
    stars: int
    created_at: datetime
    pushed_at: datetime
    html_url: str
    default_branch: str
    description: str
    archived: bool
    fork: bool
    size_kb: int


@dataclass
class WorkflowTask:
    task_id: str
    workflow_file_id: str
    repo_id: str
    name: str
    engine: str
    native_type: str
    definition: str
    start_line: int
    end_line: int


@dataclass
class Artifact:
    artifact_id: str
    workflow_file_id: str
    repo_id: str
    value: str
    artifact_type: str


def git_clone_or_update(repo: RepoCandidate, repos_dir: Path, depth: int = CLONE_DEPTH) -> Path:
    dest = repos_dir / safe_repo_dir_name(repo.full_name)
    branch = repo.default_branch

    if dest.is_dir():
        if UPDATE_EXISTING_REPOS:
            log_message(f"Updating {repo.full_name}")
            subprocess.run(
                ["git", "-C", str(dest), "fetch", "--depth", str(depth), "origin", branch],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["git", "-C", str(dest), "reset", "--hard", f"origin/{branch}"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        else:
            log_message(f"Reusing existing local clone for {repo.full_name}")
    else:
        log_message(f"Cloning {repo.full_name}")
        subprocess.run(
            [
                "git",
                "clone",
                f"--depth={depth}",
                "--branch",
                branch,
                f"https://github.com/{repo.full_name}.git",
                str(dest),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    return dest


def list_git_files(repo_path: Path) -> List[str]:
    out = subprocess.check_output(
        ["git", "-C", str(repo_path), "ls-tree", "-r", "--name-only", "HEAD"],
        stderr=subprocess.DEVNULL,
    )
    return out.decode("utf-8", errors="ignore").splitlines()


def extract_git_log(repo_path: Path, limit: int) -> List[Tuple[str, str, str, str]]:
    sep = "\x1f"
    out = subprocess.check_output(
        [
            "git",
            "-C",
            str(repo_path),
            "log",
            f"-n{limit}",
            f"--pretty=format:%H{sep}%an{sep}%ad{sep}%s",
            "--date=iso8601-strict",
        ],
        stderr=subprocess.DEVNULL,
    ).decode("utf-8", errors="replace")

    rows = []
    for line in out.splitlines():
        parts = line.split(sep, 3)
        if len(parts) != 4:
            continue
        sha, author, date_str, message = parts
        rows.append((sha, author, date_str, message))
    return rows


def build_search_query(engine: str, min_stars: int, cutoff: datetime) -> str:
    topic_part = f"topic:{engine}"
    base = [topic_part, f"stars:>={min_stars}", "fork:false", f"pushed:>={cutoff.date().isoformat()}"]
    return " ".join(base)


def to_candidate(repo) -> RepoCandidate:
    return RepoCandidate(
        repo_id=str(repo.id),
        full_name=repo.full_name,
        name=repo.name,
        owner=repo.owner.login,
        stars=repo.stargazers_count,
        created_at=repo.created_at,
        pushed_at=repo.pushed_at,
        html_url=repo.html_url,
        default_branch=repo.default_branch,
        description=repo.description or "",
        archived=repo.archived,
        fork=repo.fork,
        size_kb=repo.size,
    )


def has_engine_files(file_paths: Sequence[str], engine: str) -> bool:
    return any(matches_workflow_file(p, engine, purpose="repo") for p in file_paths)


def workflow_quality_score(repo: RepoCandidate, file_paths: Sequence[str], engine: str) -> int:
    workflow_files = sum(
        1 for p in file_paths
        if matches_workflow_file(p, engine, purpose="repo")
    )

    lower = [p.lower() for p in file_paths]

    if engine == "snakemake":
        module_bonus = sum(1 for p in lower if "rules/" in p or "workflow/" in p)
    elif engine == "nextflow":
        module_bonus = sum(1 for p in lower if "modules/" in p or "subworkflows/" in p)
    else:
        module_bonus = 0
        
    score = 0
    score += min(repo.stars, 200)
    score += workflow_files * 20
    score += min(module_bonus * 5, 30)
    score += 20 if not repo.archived else -30
    return score


def select_repositories(engine: str, max_repos: int, repos_dir: Path) -> List[RepoCandidate]:
    cutoff = now_utc() - timedelta(days=30 * MONTHS_ACTIVE)
    query = build_search_query(engine, MIN_STARS, cutoff)
    log_message(f"[{engine}] Repository search query: {query}")

    broad_candidates: List[RepoCandidate] = []
    for repo in gh.search_repositories(query, sort="stars", order="desc"):
        try:
            cand = to_candidate(repo)
            if cand.archived or cand.fork:
                continue
            broad_candidates.append(cand)
            if len(broad_candidates) >= max_repos * 4:
                break
        except GithubException as e:
            log_message(f"[{engine}] Skipping repository because of GitHub API error: {e}")
        except Exception as e:
            log_message(f"[{engine}] Skipping repository because of unexpected error: {e}")

    selected: List[Tuple[int, RepoCandidate]] = []

    for cand in broad_candidates:
        try:
            repo_path = git_clone_or_update(cand, repos_dir=repos_dir, depth=CLONE_DEPTH)
            files = list_git_files(repo_path)
            if not has_engine_files(files, engine):
                log_message(f"[{engine}] Excluded {cand.full_name}: no valid {engine} workflow files found")
                continue
            score = workflow_quality_score(cand, files, engine)
            selected.append((score, cand))
        except subprocess.CalledProcessError as e:
            log_message(f"[{engine}] Git error while validating {cand.full_name}: {e}")
        except Exception as e:
            log_message(f"[{engine}] Validation failed for {cand.full_name}: {e}")

    selected.sort(key=lambda item: (-item[0], item[1].full_name.lower()))
    final_repos = [cand for _, cand in selected[:max_repos]]
    log_message(f"[{engine}] Selected {len(final_repos)} repositories after workflow validation")
    return final_repos


def graphql_request(query: str, variables: Dict[str, object]) -> Dict[str, object]:
    response = requests.post(
        GITHUB_GRAPHQL,
        json={"query": query, "variables": variables},
        headers={"Authorization": f"bearer {TOKEN}"},
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    payload = response.json()
    if "errors" in payload:
        raise RuntimeError(f"GraphQL returned errors: {payload['errors']}")
    return payload


def fetch_issues_and_prs(repo: RepoCandidate, raw_dir: Path) -> Dict[str, object]:
    owner, name = repo.full_name.split("/", 1)

    query = """
    query($owner:String!, $name:String!, $issuesFirst:Int!, $prsFirst:Int!) {
      repository(owner:$owner, name:$name) {
        issues(first:$issuesFirst, orderBy:{field:CREATED_AT, direction:DESC}, states:[OPEN, CLOSED]) {
          nodes {
            number
            title
            state
            createdAt
            closedAt
            author { login }
          }
        }
        pullRequests(first:$prsFirst, orderBy:{field:CREATED_AT, direction:DESC}, states:[OPEN, CLOSED, MERGED]) {
          nodes {
            number
            title
            state
            createdAt
            mergedAt
            author { login }
            files(first:100) { nodes { path } }
            closingIssuesReferences(first:100) { nodes { number } }
          }
        }
      }
    }
    """

    payload = graphql_request(
        query,
        {
            "owner": owner,
            "name": name,
            "issuesFirst": ISSUES_PER_REPO,
            "prsFirst": PRS_PER_REPO,
        },
    )

    time.sleep(SLEEP_BETWEEN_GRAPHQL_CALLS)
    data = safe_get(payload, "data", "repository", default={})

    raw_path = raw_dir / "graphql" / f"{safe_repo_dir_name(repo.full_name)}.json"
    write_json(raw_path, payload)
    return data if isinstance(data, dict) else {}


def line_span(text: str, block: str) -> Tuple[int, int]:
    idx = text.find(block)
    if idx == -1:
        return -1, -1
    start_line = text.count("\n", 0, idx) + 1
    end_line = start_line + block.count("\n")
    return start_line, end_line


def extract_snakemake_io(block: str) -> Tuple[List[str], List[str], List[str]]:
    inputs: List[str] = []
    outputs: List[str] = []
    refs: List[str] = []

    m_in = SECTION_INPUT_RX.search(block)
    if m_in:
        segment = m_in.group(1)
        inputs.extend(QUOTED_STRING_RX.findall(segment))
        refs.extend(RULE_REF_RX.findall(segment))

    m_out = SECTION_OUTPUT_RX.search(block)
    if m_out:
        segment = m_out.group(1)
        outputs.extend(QUOTED_STRING_RX.findall(segment))

    return sorted(set(inputs)), sorted(set(outputs)), sorted(set(refs))


def extract_nextflow_channels(block: str) -> Tuple[List[str], List[str]]:
    inputs: List[str] = []
    outputs: List[str] = []

    m_in = SECTION_INPUT_RX.search(block)
    if m_in:
        tokens = CHANNEL_NAME_RX.findall(m_in.group(1))
        inputs.extend(tokens)

    m_out = SECTION_OUTPUT_RX.search(block)
    if m_out:
        tokens = CHANNEL_NAME_RX.findall(m_out.group(1))
        outputs.extend(tokens)

    stop = {
        "path", "val", "tuple", "env", "stdin", "stdout", "emit", "into", "from",
        "file", "script", "shell", "exec", "when", "each", "true", "false", "input", "output"
    }
    inputs = [t for t in inputs if t not in stop and len(t) > 1]
    outputs = [t for t in outputs if t not in stop and len(t) > 1]
    return sorted(set(inputs)), sorted(set(outputs))


def open_csv_writers(data_dir: Path) -> Tuple[Dict[str, csv.writer], Dict[str, object]]:
    specs = {
        "repos": ("repos.csv", [
            "repo_id", "name", "full_name", "owner", "stars", "created_at", "pushed_at",
            "html_url", "default_branch", "description", "archived", "size_kb", "engine"
        ]),
        "files": ("files.csv", ["file_id", "repo_id", "path", "extension", "is_workflow_file"]),
        "commits": ("commits.csv", ["repo_id", "sha", "author", "date", "message"]),
        "issues": ("issues.csv", ["repo_id", "issue_id", "author", "title", "state", "created_at", "closed_at"]),
        "prs": ("prs.csv", ["repo_id", "pr_id", "author", "title", "state", "created_at", "merged_at"]),
        "pr_files": ("pr_files.csv", ["repo_id", "pr_id", "file_id"]),
        "pr_closes": ("pr_closes.csv", ["repo_id", "pr_id", "issue_id"]),
        "workflow_files": ("workflow_files.csv", ["workflow_file_id", "repo_id", "file_id", "path", "engine", "file_type"]),
        "tasks": ("tasks.csv", [
            "task_id", "workflow_file_id", "repo_id", "name", "engine", "native_type",
            "definition", "start_line", "end_line"
        ]),
        "artifacts": ("artifacts.csv", ["artifact_id", "workflow_file_id", "repo_id", "value", "artifact_type"]),
        "task_reads": ("task_reads.csv", ["task_id", "artifact_id"]),
        "task_writes": ("task_writes.csv", ["task_id", "artifact_id"]),
        "task_dependencies": ("task_dependencies.csv", ["src_task_id", "dst_task_id", "dependency_type", "evidence"]),
    }

    writers: Dict[str, csv.writer] = {}
    handles: Dict[str, object] = {}

    for key, (filename, header) in specs.items():
        path = data_dir / filename
        fh = path.open("w", newline="", encoding="utf-8")
        writer = csv.writer(fh)
        writer.writerow(header)
        writers[key] = writer
        handles[key] = fh

    return writers, handles


def run_for_engine(engine: str) -> None:
    repos_dir, data_dir, raw_dir = ensure_dirs(engine)
    writers, handles = open_csv_writers(data_dir)

    selected = select_repositories(engine, MAX_REPOS, repos_dir)

    artifact_seen: set[str] = set()
    file_map: Dict[Tuple[str, str], str] = {}
    workflow_file_map: Dict[Tuple[str, str], str] = {}

    for repo in selected:
        repo_path = repos_dir / safe_repo_dir_name(repo.full_name)

        writers["repos"].writerow([
            repo.repo_id,
            repo.name,
            repo.full_name,
            repo.owner,
            repo.stars,
            repo.created_at.isoformat(),
            repo.pushed_at.isoformat(),
            repo.html_url,
            repo.default_branch,
            repo.description.replace("\n", " "),
            str(repo.archived).lower(),
            repo.size_kb,
            engine,
        ])

        try:
            for sha, author, date_str, message in extract_git_log(repo_path, COMMITS_PER_REPO):
                writers["commits"].writerow([repo.repo_id, sha, author, date_str, message])
        except Exception as e:
            log_message(f"[{engine}] Commit extraction failed for {repo.full_name}: {e}")

        try:
            files = list_git_files(repo_path)
        except Exception as e:
            log_message(f"[{engine}] File listing failed for {repo.full_name}: {e}")
            continue

        task_name_to_id: Dict[str, str] = {}
        snakemake_refs_buffer: List[Tuple[str, str]] = []
        nextflow_channel_reads: Dict[str, List[str]] = {}
        nextflow_channel_writes: Dict[str, List[str]] = {}

        for rel_path in files:
            file_id = stable_id(repo.repo_id, rel_path)
            ext = Path(rel_path).suffix.lower().lstrip(".")
            is_workflow_file = matches_workflow_file(rel_path, engine, purpose="repo")
            file_map[(repo.repo_id, rel_path)] = file_id

            writers["files"].writerow([
                file_id,
                repo.repo_id,
                rel_path,
                ext,
                str(is_workflow_file).lower(),
            ])

            full_path = repo_path / rel_path
            if not full_path.is_file():
                continue

            is_snakemake = matches_workflow_file(rel_path, "snakemake", purpose="task")
            is_nextflow = matches_workflow_file(rel_path, "nextflow", purpose="task")

            if not (is_snakemake or is_nextflow):
                continue

            workflow_file_id = stable_id(repo.repo_id, "workflow", rel_path)
            workflow_file_map[(repo.repo_id, rel_path)] = workflow_file_id
            file_type = "workflow_source"
            workflow_engine = "snakemake" if is_snakemake else "nextflow"

            writers["workflow_files"].writerow([
                workflow_file_id, repo.repo_id, file_id, rel_path, workflow_engine, file_type
            ])

            text = read_text_file(full_path)
            raw_text_path = raw_dir / "workflow_text" / safe_repo_dir_name(repo.full_name) / rel_path.replace("/", "__")
            raw_text_path.parent.mkdir(parents=True, exist_ok=True)
            raw_text_path.write_text(text, encoding="utf-8")

            if is_snakemake:
                matches: List[Tuple[str, str, str]] = []
                matches.extend(("rule", name, body) for name, body in RULE_RX.findall(text))
                matches.extend(("checkpoint", name, body) for name, body in CHECKPOINT_RX.findall(text))

                for native_type, name, body in matches:
                    full_block = f"{native_type} {name}:" + body
                    start_line, end_line = line_span(text, full_block)
                    task_id = stable_id(repo.repo_id, rel_path, native_type, name)
                    task_name_to_id[name] = task_id

                    clean_definition = sanitize_for_csv(body)
                    writers["tasks"].writerow([
                        task_id,
                        workflow_file_id,
                        repo.repo_id,
                        name,
                        "snakemake",
                        native_type,
                        clean_definition,
                        start_line,
                        end_line,
                    ])

                    ins, outs, refs = extract_snakemake_io(body)
                    for value in ins:
                        artifact_id = stable_id(repo.repo_id, workflow_file_id, "artifact", value)
                        if artifact_id not in artifact_seen:
                            artifact_seen.add(artifact_id)
                            writers["artifacts"].writerow([artifact_id, workflow_file_id, repo.repo_id, value, "file_pattern"])
                        writers["task_reads"].writerow([task_id, artifact_id])

                    for value in outs:
                        artifact_id = stable_id(repo.repo_id, workflow_file_id, "artifact", value)
                        if artifact_id not in artifact_seen:
                            artifact_seen.add(artifact_id)
                            writers["artifacts"].writerow([artifact_id, workflow_file_id, repo.repo_id, value, "file_pattern"])
                        writers["task_writes"].writerow([task_id, artifact_id])

                    for ref_name in refs:
                        snakemake_refs_buffer.append((task_id, ref_name))

            if is_nextflow:
                for name, body in PROCESS_RX.findall(text):
                    full_block = f"process {name} {{" + body
                    start_line, end_line = line_span(text, full_block)
                    task_id = stable_id(repo.repo_id, rel_path, "process", name)
                    task_name_to_id[name] = task_id

                    clean_definition = sanitize_for_csv(body)
                    writers["tasks"].writerow([
                        task_id,
                        workflow_file_id,
                        repo.repo_id,
                        name,
                        "nextflow",
                        "process",
                        clean_definition,
                        start_line,
                        end_line,
                    ])

                    ins, outs = extract_nextflow_channels(body)
                    nextflow_channel_reads[task_id] = ins
                    nextflow_channel_writes[task_id] = outs

                    for value in ins:
                        artifact_id = stable_id(repo.repo_id, workflow_file_id, "channel", value)
                        if artifact_id not in artifact_seen:
                            artifact_seen.add(artifact_id)
                            writers["artifacts"].writerow([artifact_id, workflow_file_id, repo.repo_id, value, "channel"])
                        writers["task_reads"].writerow([task_id, artifact_id])

                    for value in outs:
                        artifact_id = stable_id(repo.repo_id, workflow_file_id, "channel", value)
                        if artifact_id not in artifact_seen:
                            artifact_seen.add(artifact_id)
                            writers["artifacts"].writerow([artifact_id, workflow_file_id, repo.repo_id, value, "channel"])
                        writers["task_writes"].writerow([task_id, artifact_id])

                # channel_producers: Dict[str, List[str]] = {}
                # channel_consumers: Dict[str, List[str]] = {}
                # for task_id, channels in nextflow_channel_writes.items():
                #     for ch in channels:
                #         channel_producers.setdefault(ch, []).append(task_id)
                # for task_id, channels in nextflow_channel_reads.items():
                #     for ch in channels:
                #         channel_consumers.setdefault(ch, []).append(task_id)

                # for ch, producers in channel_producers.items():
                #     for producer in producers:
                #         for consumer in channel_consumers.get(ch, []):
                #             if producer != consumer:
                #                 writers["task_dependencies"].writerow([
                #                     consumer,
                #                     producer,
                #                     "channel_flow",
                #                     ch,
                #                 ])

        for src_task_id, ref_name in snakemake_refs_buffer:
            target_id = task_name_to_id.get(ref_name)
            if target_id and target_id != src_task_id:
                writers["task_dependencies"].writerow([
                    src_task_id,
                    target_id,
                    "explicit_rule_reference",
                    ref_name,
                ])

        try:
            data = fetch_issues_and_prs(repo, raw_dir)
        except Exception as e:
            log_message(f"[{engine}] Issue/PR extraction failed for {repo.full_name}: {e}")
            continue

        for issue in safe_nodes(data, "issues", "nodes"):
            login = safe_get(issue, "author", "login", default="")
            writers["issues"].writerow([
                repo.repo_id,
                issue.get("number", ""),
                login,
                (issue.get("title") or "").replace("\n", " "),
                issue.get("state", ""),
                issue.get("createdAt", ""),
                issue.get("closedAt") or "",
            ])

        for pr in safe_nodes(data, "pullRequests", "nodes"):
            login = safe_get(pr, "author", "login", default="")
            pr_id = pr.get("number", "")
            writers["prs"].writerow([
                repo.repo_id,
                pr_id,
                login,
                (pr.get("title") or "").replace("\n", " "),
                pr.get("state", ""),
                pr.get("createdAt", ""),
                pr.get("mergedAt") or "",
            ])

            pr_file_nodes = safe_nodes(pr, "files", "nodes")
            if not pr_file_nodes:
                log_message(f"[{engine}] Skipping PR file links for {repo.full_name} PR #{pr_id}: no file list returned")

            for changed in pr_file_nodes:
                changed_path = changed.get("path")
                if not changed_path:
                    continue
                file_id = file_map.get((repo.repo_id, changed_path))
                if file_id:
                    writers["pr_files"].writerow([repo.repo_id, pr_id, file_id])

            for closed in safe_nodes(pr, "closingIssuesReferences", "nodes"):
                issue_number = closed.get("number")
                if issue_number is not None:
                    writers["pr_closes"].writerow([repo.repo_id, pr_id, issue_number])

    for fh in handles.values():
        fh.close()

    log_message(f"[{engine}] Extraction complete. CSV files written to: {data_dir.resolve()}")


def main() -> None:
    engines = ["snakemake", "nextflow"] if ENGINE == "both" else [ENGINE]

    for engine in engines:
        log_message(f"=== Starting extraction for engine: {engine} ===")
        run_for_engine(engine)
        log_message(f"=== Finished extraction for engine: {engine} ===")


if __name__ == "__main__":
    main()