# CI/CD Pipeline

This document provides an overview of the Continuous Integration and Continuous Deployment pipeline for the ASD Treatment Database.

## Overview

The project uses GitHub Actions to automate testing. The workflow is defined in `.github/workflows/ci-cd.yml`.

## Workflow Triggers

- Push to `main` branch
- Pull requests to `main` branch
- Manual workflow dispatch via GitHub Actions UI

## Test Job

The automated test job performs the following:

### Setup
- Python 3.10
- Node.js 18
- Dependency caching for faster builds

### Testing
- Installs minimal dependencies from `backend/requirements-test.txt`
- Runs backend pytest with `DISABLE_MODEL_LOADING=1` to skip heavy ML model loading
- Excludes end-to-end and performance tests for speed

### Purpose of DISABLE_MODEL_LOADING

This environment variable prevents loading MedBERT and sentence transformer models during tests, significantly reducing:
- Test execution time
- Memory usage
- Dependency download size

Mock objects are used instead for testing API logic.

## Deploy Job

**Current Status:** Deployment is currently disabled (`if: false`) in the workflow.

The deployment strategy is being updated to use a non-SSH approach. Documentation will be updated once the new deployment method is implemented.

## Manual Deployment

To deploy manually to production:

```bash
cd /path/to/project
git fetch origin
git reset --hard origin/main
./deploy-production.sh
```

The `deploy-production.sh` script handles:
- Environment setup
- Frontend build
- Service restarts
- Nginx reload
- Health checks

## Workflow Concurrency

The workflow prevents multiple simultaneous runs:

```yaml
concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true
```

## Best Practices

### Before Merging
- Run tests locally: `./run_all_tests.sh`
- Verify tests pass in GitHub Actions
- Get code review approval

### Branch Protection
Configure on `main` branch:
- Require pull request reviews
- Require status checks to pass
- Require conversation resolution

### Version Tagging
Tag releases for tracking:

```bash
git tag -a v1.0.0 -m "Release version 1.0.0"
git push origin v1.0.0
```
