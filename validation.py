#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path

DATA_DIR = Path("data")

EXPECTED_FILES = [
    "repos.csv",
    "files.csv",
    "workflow_files.csv",
    "tasks.csv",
    "artifacts.csv",
    "issues.csv",
    "prs.csv",
    "pr_files.csv",
]

def count_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as f:
        return max(sum(1 for _ in f) - 1, 0)

def read_csv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        yield from csv.DictReader(f)

def describe_counter(name: str, counter: Counter):
    values = list(counter.values())
    if not values:
        print(f"{name}: no data")
        return
    print(f"{name}:")
    print(f"  repos with data: {len(values)}")
    print(f"  min: {min(values)}")
    print(f"  max: {max(values)}")
    print(f"  avg: {sum(values)/len(values):.2f}")
    print("  top 10:")
    for repo_id, n in counter.most_common(10):
        print(f"    {repo_id}: {n}")

def load_repo_names(engine_dir: Path) -> dict[str, str]:
    repos_path = engine_dir / "repos.csv"
    if not repos_path.exists():
        return {}
    repos = list(read_csv(repos_path))
    return {r["repo_id"]: r["full_name"] for r in repos}

def named(counter: Counter, repo_names: dict[str, str]) -> Counter:
    out = Counter()
    for repo_id, n in counter.items():
        out[repo_names.get(repo_id, repo_id)] = n
    return out

def count_by_repo(path: Path) -> Counter:
    counter = Counter()
    if not path.exists():
        return counter
    for row in read_csv(path):
        repo_id = row.get("repo_id")
        if repo_id:
            counter[repo_id] += 1
    return counter

def validate_engine(engine_dir: Path):
    print(f"\n{'=' * 60}")
    print(f"ENGINE: {engine_dir.name}")
    print(f"{'=' * 60}")

    print("=== Row counts ===")
    for filename in EXPECTED_FILES:
        path = engine_dir / filename
        if path.exists():
            print(f"{filename}: {count_rows(path)}")
        else:
            print(f"{filename}: MISSING")

    repo_names = load_repo_names(engine_dir)

    files_per_repo = count_by_repo(engine_dir / "files.csv")
    workflow_files_per_repo = count_by_repo(engine_dir / "workflow_files.csv")
    tasks_per_repo = count_by_repo(engine_dir / "tasks.csv")
    artifacts_per_repo = count_by_repo(engine_dir / "artifacts.csv")
    issues_per_repo = count_by_repo(engine_dir / "issues.csv")
    prs_per_repo = count_by_repo(engine_dir / "prs.csv")
    pr_files_per_repo = count_by_repo(engine_dir / "pr_files.csv")

    print("\n=== Per-repo distributions ===")
    describe_counter("files per repo", named(files_per_repo, repo_names))
    describe_counter("workflow files per repo", named(workflow_files_per_repo, repo_names))
    describe_counter("tasks per repo", named(tasks_per_repo, repo_names))
    describe_counter("artifacts per repo", named(artifacts_per_repo, repo_names))
    describe_counter("issues per repo", named(issues_per_repo, repo_names))
    describe_counter("PRs per repo", named(prs_per_repo, repo_names))
    describe_counter("PR-file links per repo", named(pr_files_per_repo, repo_names))

def main():
    engine_dirs = sorted([p for p in DATA_DIR.iterdir() if p.is_dir()])

    if not engine_dirs:
        print(f"No engine subdirectories found in {DATA_DIR}")
        return

    for engine_dir in engine_dirs:
        validate_engine(engine_dir)

if __name__ == "__main__":
    main()