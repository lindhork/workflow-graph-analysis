# Graph-Based Exploration of Scientific Workflow Projects on GitHub

This repository contains the implementation developed for the study project  
**"Graph-Based Exploration of Scientific Workflow Projects on GitHub"**.

The project investigates how scientific workflow repositories can be analyzed
using a graph-based representation that integrates workflow structure and
software development metadata.

## Project Overview

Scientific Workflow Management Systems such as **Snakemake** and **Nextflow**
allow researchers to describe computational pipelines as directed acyclic
graphs (DAGs). At the same time, many workflow projects are hosted on
GitHub, where repositories also contain development metadata such as
commits, issues, and pull requests.

This project combines both perspectives by constructing a **property graph**
that integrates:

- workflow structure (tasks, artifacts, dependencies)
- repository structure (files and workflow definitions)
- development metadata (commits, issues, pull requests)

The resulting graph representation can be explored using **Neo4j** and
Cypher queries.

## Repository Contents

The repository includes the following components:

- **GitHub extraction pipeline**  
  Python scripts used to collect repository metadata and workflow files
  from GitHub.

- **CSV export**  
  Structured data export used for graph database import.

- **Neo4j import scripts**  
  Cypher scripts used to import the dataset into a Neo4j graph database.

## Data Model

The graph model integrates workflow entities and repository metadata using
the following main node types:

- Repository
- File
- WorkflowFile
- Task
- Artifact
- Commit
- Issue
- PullRequest
- Contributor

Relationships capture workflow dependencies, artifact usage, repository
structure, and development activity.

## Reproducibility

To reproduce the analysis:

1. Run the extraction pipeline to collect repository data.
3. Import the dataset into Neo4j using the provided import scripts.
4. Execute Cypher queries to explore workflow structures and repository activity.

## Technologies Used

- Python
- GitHub API
- Neo4j Graph Database
- Cypher Query Language

