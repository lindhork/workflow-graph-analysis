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

- **workflow_github_extractor.py**  
  Python script used to collect repository metadata and workflow files
  from GitHub and export the extracted data as CSV files.

- **validation.py**  
  Script used to inspect the generated CSV files and summarize dataset
  statistics per workflow engine and repository.

- **import.cypher**  
  Cypher script used to import the generated CSV files into the Neo4j
  graph database.

- **requirements.txt**  
  Python dependencies required to run the extraction pipeline.

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

## Requirements

- Python >= 3.10
- GitHub API access
- Neo4j Graph Database

## Setup

Install the required Python dependencies:

    pip install -r requirements.txt

### GitHub API Token

The extraction pipeline requires a GitHub personal access token to access
the GitHub API and avoid strict rate limits.

Set the token as an environment variable:

    export GITHUB_TOKEN=<your_github_token>

## Usage

Run the extraction pipeline to collect repository data:

    python workflow_github_extractor.py

The script generates CSV files containing repositories, workflow tasks,
artifacts, files, commits, issues, and pull requests.

## Importing Data into Neo4j

After generating the CSV files, import them into Neo4j using the provided
Cypher script.

Open the Neo4j Browser and execute:

    :source import.cypher

This will create the graph model and load all nodes and relationships
into the database.

## Validation

The repository also includes a validation script that provides a basic
sanity check for the generated CSV exports.

Run:

    python validation.py

The script expects engine-specific subdirectories inside `data/` and reports:

- whether the expected CSV files are present
- row counts for each CSV file
- per-repository distributions for files, workflow files, tasks, artifacts,
  issues, pull requests, and pull request–file links

This helps to inspect the extracted dataset before importing it into Neo4j.

## Technologies Used

- Python
- GitHub API
- PyGithub
- Neo4j Graph Database
- Cypher Query Language
