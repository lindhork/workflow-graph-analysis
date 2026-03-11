
CREATE CONSTRAINT repo_id_unique IF NOT EXISTS
FOR (n:Repository) REQUIRE n.repo_id IS UNIQUE;

CREATE CONSTRAINT file_id_unique IF NOT EXISTS
FOR (n:File) REQUIRE n.file_id IS UNIQUE;

CREATE CONSTRAINT workflow_file_id_unique IF NOT EXISTS
FOR (n:WorkflowFile) REQUIRE n.workflow_file_id IS UNIQUE;

CREATE CONSTRAINT task_id_unique IF NOT EXISTS
FOR (n:Task) REQUIRE n.task_id IS UNIQUE;

CREATE CONSTRAINT artifact_id_unique IF NOT EXISTS
FOR (n:Artifact) REQUIRE n.artifact_id IS UNIQUE;

CREATE CONSTRAINT commit_sha_unique IF NOT EXISTS
FOR (n:Commit) REQUIRE n.sha IS UNIQUE;

CREATE CONSTRAINT issue_key_unique IF NOT EXISTS
FOR (n:Issue) REQUIRE n.issue_key IS UNIQUE;

CREATE CONSTRAINT pr_key_unique IF NOT EXISTS
FOR (n:PullRequest) REQUIRE n.pr_key IS UNIQUE;

CREATE CONSTRAINT contributor_login_unique IF NOT EXISTS
FOR (n:Contributor) REQUIRE n.login IS UNIQUE;


LOAD CSV WITH HEADERS FROM 'file:///repos.csv' AS row
MERGE (r:Repository {repo_id: row.repo_id})
SET r.name = row.name,
    r.full_name = row.full_name,
    r.owner = row.owner,
    r.stars = toInteger(row.stars),
    r.created_at = datetime(row.created_at),
    r.pushed_at = datetime(row.pushed_at),
    r.html_url = row.html_url,
    r.default_branch = row.default_branch,
    r.description = row.description,
    r.archived = row.archived = 'true',
    r.size_kb = toInteger(row.size_kb),
    r.engine = row.engine;


LOAD CSV WITH HEADERS FROM 'file:///files.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
MERGE (f:File {file_id: row.file_id})
SET f.path = row.path,
    f.extension = row.extension,
    f.is_workflow_file = row.is_workflow_file = 'true'
MERGE (r)-[:CONTAINS_FILE]->(f);


LOAD CSV WITH HEADERS FROM 'file:///workflow_files.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
MATCH (f:File {file_id: row.file_id})
MERGE (wf:WorkflowFile {workflow_file_id: row.workflow_file_id})
SET wf.path = row.path,
    wf.engine = row.engine,
    wf.file_type = row.file_type
MERGE (r)-[:CONTAINS_WORKFLOW_FILE]->(wf)
MERGE (wf)-[:REPRESENTS_FILE]->(f);


LOAD CSV WITH HEADERS FROM 'file:///tasks.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
MATCH (wf:WorkflowFile {workflow_file_id: row.workflow_file_id})
MERGE (t:Task {task_id: row.task_id})
SET t.name = row.name,
    t.engine = row.engine,
    t.native_type = row.native_type,
    t.definition = row.definition,
    t.start_line = CASE WHEN row.start_line = '' THEN null ELSE toInteger(row.start_line) END,
    t.end_line = CASE WHEN row.end_line = '' THEN null ELSE toInteger(row.end_line) END
MERGE (wf)-[:DECLARES_TASK]->(t)
MERGE (r)-[:CONTAINS_TASK]->(t);


LOAD CSV WITH HEADERS FROM 'file:///artifacts.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
MATCH (wf:WorkflowFile {workflow_file_id: row.workflow_file_id})
MERGE (a:Artifact {artifact_id: row.artifact_id})
SET a.value = row.value,
    a.artifact_type = row.artifact_type
MERGE (wf)-[:DECLARES_ARTIFACT]->(a)
MERGE (r)-[:USES_ARTIFACT]->(a);


LOAD CSV WITH HEADERS FROM 'file:///task_reads.csv' AS row
MATCH (t:Task {task_id: row.task_id})
MATCH (a:Artifact {artifact_id: row.artifact_id})
MERGE (t)-[:READS]->(a);


LOAD CSV WITH HEADERS FROM 'file:///task_writes.csv' AS row
MATCH (t:Task {task_id: row.task_id})
MATCH (a:Artifact {artifact_id: row.artifact_id})
MERGE (t)-[:WRITES]->(a);


LOAD CSV WITH HEADERS FROM 'file:///task_dependencies.csv' AS row
MATCH (src:Task {task_id: row.src_task_id})
MATCH (dst:Task {task_id: row.dst_task_id})
MERGE (src)-[d:DEPENDS_ON]->(dst)
SET d.dependency_type = row.dependency_type,
    d.evidence = row.evidence;


LOAD CSV WITH HEADERS FROM 'file:///commits.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
MERGE (c:Commit {sha: row.sha})
SET c.author = row.author,
    c.date = datetime(row.date),
    c.message = row.message
MERGE (r)-[:HAS_COMMIT]->(c)
FOREACH (_ IN CASE WHEN row.author IS NULL OR trim(row.author) = '' THEN [] ELSE [1] END |
  MERGE (u:Contributor {login: row.author})
  MERGE (u)-[:AUTHORED_COMMIT]->(c)
);


LOAD CSV WITH HEADERS FROM 'file:///issues.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
WITH r, row, row.repo_id + '::' + row.issue_id AS issue_key
MERGE (i:Issue {issue_key: issue_key})
SET i.issue_id = row.issue_id,
    i.author = row.author,
    i.title = row.title,
    i.state = row.state,
    i.created_at = datetime(row.created_at),
    i.closed_at = CASE WHEN row.closed_at = '' THEN null ELSE datetime(row.closed_at) END
MERGE (r)-[:HAS_ISSUE]->(i);


LOAD CSV WITH HEADERS FROM 'file:///prs.csv' AS row
MATCH (r:Repository {repo_id: row.repo_id})
WITH r, row, row.repo_id + '::' + row.pr_id AS pr_key
MERGE (p:PullRequest {pr_key: pr_key})
SET p.pr_id = row.pr_id,
    p.author = row.author,
    p.title = row.title,
    p.state = row.state,
    p.created_at = datetime(row.created_at),
    p.merged_at = CASE WHEN row.merged_at = '' THEN null ELSE datetime(row.merged_at) END
MERGE (r)-[:HAS_PR]->(p);


LOAD CSV WITH HEADERS FROM 'file:///pr_files.csv' AS row
MATCH (p:PullRequest {pr_key: row.repo_id + '::' + row.pr_id})
MATCH (f:File {file_id: row.file_id})
MERGE (p)-[:MODIFIES]->(f);


LOAD CSV WITH HEADERS FROM 'file:///pr_closes.csv' AS row
MATCH (p:PullRequest {pr_key: row.repo_id + '::' + row.pr_id})
MATCH (i:Issue {issue_key: row.repo_id + '::' + row.issue_id})
MERGE (p)-[:CLOSES]->(i);
