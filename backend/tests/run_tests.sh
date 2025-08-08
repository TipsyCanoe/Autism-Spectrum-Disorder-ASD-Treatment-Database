#!/bin/bash
# run_tests.sh - Script to run tests with coverage reporting

# Define color codes for prettier output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Help message
show_help() {
    echo -e "${BLUE}Usage:${NC} $0 [options]"
    echo -e ""
    echo -e "${BLUE}Options:${NC}"
    echo -e "  -h, --help         Show this help message"
    echo -e "  -m, --marker       Run tests with specific marker (unit, api, integration, e2e, performance)"
    echo -e "  -f, --file         Run tests from specific file"
    echo -e "  -k, --keyword      Run tests matching keyword expression"
    echo -e "  --no-coverage      Run tests without coverage reporting"
    echo -e "  -p, --parallel     Run tests in parallel (faster for large test suites)"
    echo -e "  --html-report      Generate HTML coverage report"
    echo -e ""
    echo -e "${BLUE}Examples:${NC}"
    echo -e "  $0                 Run all tests with coverage"
    echo -e "  $0 -m unit         Run only unit tests"
    echo -e "  $0 -f test_api.py  Run only tests in test_api.py"
    echo -e "  $0 -k cache        Run tests with 'cache' in the name"
    echo -e "  $0 --no-coverage   Run all tests without coverage"
    echo -e "  $0 -p              Run tests in parallel"
    echo -e "  $0 --html-report   Generate HTML coverage report"
}

# Default values
marker=""
file=""
keyword=""
coverage=true
parallel=false
html_report=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -m|--marker)
            marker="$2"
            shift 2
            ;;
        -f|--file)
            file="$2"
            shift 2
            ;;
        -k|--keyword)
            keyword="$2"
            shift 2
            ;;
        --no-coverage)
            coverage=false
            shift
            ;;
        -p|--parallel)
            parallel=true
            shift
            ;;
        --html-report)
            html_report=true
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

echo -e "${YELLOW}Installing required packages...${NC}"

# Install pytest and coverage packages if not installed
pip install pytest pytest-cov || pip3 install pytest pytest-cov

# Locate and activate the project's top-level virtual environment
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_PATH")")"
VENV_ACTIVATE="$PROJECT_ROOT/venv/bin/activate"
# Activate venv if present, otherwise warn and continue using system python
if [ -f "$VENV_ACTIVATE" ]; then
    source "$VENV_ACTIVATE"
else
    echo -e "${YELLOW}WARNING: venv not found at $VENV_ACTIVATE, proceeding with system python${NC}"
fi

# Locate and activate the project's virtual environment
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_PATH")")"
VENV_ACTIVATE="$PROJECT_ROOT/venv/bin/activate"
if [ -f "$VENV_ACTIVATE" ]; then
    source "$VENV_ACTIVATE"
else
    echo -e "${RED}ERROR: venv not found at $VENV_ACTIVATE${NC}"
    exit 1
fi

## Change to the backend directory to scope pytest to only backend/tests
BACKEND_DIR="$(dirname "$(dirname "${SCRIPT_PATH}")")/backend"
cd "$BACKEND_DIR" || { echo "Failed to cd to $BACKEND_DIR"; exit 1; }

# Determine which Python command to use (python or python3)
if command -v python3 &>/dev/null; then
    python_cmd="python3"
elif command -v python &>/dev/null; then
    python_cmd="python"
else
    echo -e "${RED}Error: Neither python nor python3 command found.${NC}"
    exit 1
fi

# Build the pytest command to run only backend/tests
pytest_cmd="$python_cmd -m pytest tests"

# Add marker if specified
if [[ -n "$marker" ]]; then
    pytest_cmd="$pytest_cmd -m $marker"
fi

# Add file if specified
if [[ -n "$file" ]]; then
    pytest_cmd="$pytest_cmd tests/$file"
fi

# Add keyword if specified
if [[ -n "$keyword" ]]; then
    pytest_cmd="$pytest_cmd -k $keyword"
fi

# Add parallel option if enabled
if [[ "$parallel" = true ]]; then
    echo -e "${YELLOW}Running tests in parallel...${NC}"
    pytest_cmd="$pytest_cmd -xvs"
fi

# Add coverage if enabled
if [[ "$coverage" = true ]]; then
    echo -e "${YELLOW}Running tests with coverage...${NC}"
    cov_options="--cov=.. --cov-report=term-missing:skip-covered"
    
    if [[ "$html_report" = true ]]; then
        cov_options="$cov_options --cov-report=html:coverage_report"
        echo -e "${YELLOW}HTML coverage report will be generated in coverage_report directory${NC}"
    fi
    
    pytest_cmd="$pytest_cmd $cov_options -v"
else
    echo -e "${YELLOW}Running tests without coverage...${NC}"
    pytest_cmd="$pytest_cmd -v"
fi

# Run the tests
echo -e "${BLUE}Running command:${NC} $pytest_cmd"
eval $pytest_cmd

# Show result summary
exit_code=$?
if [[ $exit_code -eq 0 ]]; then
    echo -e "${GREEN}All tests passed!${NC}"
else
    echo -e "${RED}Some tests failed.${NC}"
fi

# Return the exit code from pytest
exit $exit_code

# Get the exit status
status=$?

if [ $status -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
else
    echo -e "${RED}Some tests failed.${NC}"
fi

if [[ "$coverage" = true ]]; then
    echo -e "${YELLOW}Coverage report generated in coverage_report directory${NC}"
    echo -e "${YELLOW}Open coverage_report/index.html in a browser to view detailed results${NC}"
fi

# Run by marker examples
echo -e "\n${YELLOW}Example commands:${NC}"
echo -e "./run_tests.sh -m unit       ${GREEN}# Run only unit tests${NC}"
echo -e "./run_tests.sh -m api        ${GREEN}# Run only API tests${NC}"
echo -e "./run_tests.sh -m integration ${GREEN}# Run only integration tests${NC}"
echo -e "./run_tests.sh -m e2e        ${GREEN}# Run only end-to-end tests${NC}"
echo -e "./run_tests.sh -m performance ${GREEN}# Run only performance tests${NC}"
echo -e "./run_tests.sh -f test_api.py ${GREEN}# Run only tests in test_api.py${NC}"
echo -e "./run_tests.sh -k cache      ${GREEN}# Run tests with 'cache' in the name${NC}"

exit $status
