# Backend API Testing

This directory contains tests for the Flask backend API. The tests are organized into different categories using pytest markers.   ``pip install -U pytest`` for functionality

## Test Structure

- `conftest.py`: Contains pytest fixtures used across test files
- `pytest.ini`: Configuration file for pytest with custom markers
- `test_helpers.py`: Unit tests for helper functions
- `test_api.py`: Tests for API endpoints
- `test_cache.py`: Tests for the SimpleCache implementation
- `test_integration.py`: Integration tests that verify different components work together
- `test_e2e.py`: End-to-end tests that simulate complete user workflows
- `test_performance.py`: Tests that measure and validate performance characteristics
- `run_tests.sh`: Script to run tests with coverage reporting

## Test Categories

Tests are organized into the following categories using pytest markers:

- **unit**: Unit tests that don't require a database connection
- **api**: Tests that make HTTP requests to the API endpoints
- **integration**: Tests that verify different components work together
- **e2e**: End-to-end tests that simulate complete user workflows
- **performance**: Tests that measure and validate performance characteristics

## Running Tests

### Using the test runner script

The easiest way to run tests is using the `run_tests.sh` script, which supports several options:

```bash
# Run all tests with coverage reporting
./run_tests.sh

# Run only unit tests
./run_tests.sh -m unit

# Run tests from a specific file
./run_tests.sh -f test_api.py

# Run tests with a specific keyword in the name
./run_tests.sh -k cache

# Run tests without coverage reporting
./run_tests.sh --no-coverage

# Show help and all options
./run_tests.sh --help
```

### Run tests with verbose output:

```bash
# Using python
python -m pytest -v

# Using python3
python3 -m pytest -v
```

### Run a specific test file:

```bash
# Using python
python -m pytest tests/test_api.py

# Using python3
python3 -m pytest tests/test_api.py
```

### Run a specific test function:

```bash
# Using python
python -m pytest tests/test_api.py::test_get_filters_endpoint

# Using python3
python3 -m pytest tests/test_api.py::test_get_filters_endpoint
```

## Mock Objects

The tests use several mock objects to avoid actual database connections and external dependencies:

- `mock_db_connection`: Mocks the PostgreSQL database connection
- `mock_sentence_model`: Mocks the sentence transformer model for embeddings

## Adding New Tests

When adding new tests:

1. Use appropriate pytest markers to categorize tests
2. Use fixtures from `conftest.py` to set up test environment
3. Mock any external dependencies
4. Follow the existing test patterns for consistency
