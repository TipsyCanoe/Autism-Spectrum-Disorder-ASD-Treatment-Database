#!/bin/bash
# filepath: /home/coleoliva/senior-proj/run_all_tests.sh

# Define color codes for prettier output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Help message
show_help() {
    echo -e "${BLUE}Usage:${NC} $0 [options]"
    echo -e ""
    echo -e "${BLUE}Options:${NC}"
    echo -e "  -h, --help         Show this help message"
    echo -e "  -b, --backend      Run only backend tests"
    echo -e "  -f, --frontend     Run only frontend tests"
    echo -e "  --coverage         Generate coverage reports"
    echo -e "  --backend-marker   Run backend tests with specific marker (unit, api, integration, e2e, performance)"
    echo -e ""
    echo -e "${BLUE}Examples:${NC}"
    echo -e "  $0                 Run all tests (backend and frontend)"
    echo -e "  $0 -b              Run only backend tests"
    echo -e "  $0 -f              Run only frontend tests"
    echo -e "  $0 --coverage      Run all tests with coverage reporting"
    echo -e "  $0 -b --backend-marker unit  Run only backend unit tests"
}

# Default values
run_backend=true
run_frontend=true
coverage=false
backend_marker=""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -b|--backend)
            run_backend=true
            run_frontend=false
            shift
            ;;
        -f|--frontend)
            run_frontend=true
            run_backend=false
            shift
            ;;
        --coverage)
            coverage=true
            shift
            ;;
        --backend-marker)
            backend_marker="$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# Determine which Python command to use (python or python3)
if command -v python3 &>/dev/null; then
    python_cmd="python3"
elif command -v python &>/dev/null; then
    python_cmd="python"
else
    echo -e "${RED}Error: Neither python nor python3 command found.${NC}"
    exit 1
fi

# Function to run backend tests
run_backend_tests() {
    echo -e "\n${BOLD}${BLUE}=== Running Backend Tests ===${NC}\n"
    
    # Check if backend test directory exists
    if [ ! -d "${SCRIPT_DIR}/services/api/tests" ]; then
        echo -e "${RED}Backend test directory not found at ${SCRIPT_DIR}/services/api/tests${NC}"
        return 1
    fi
    
    # Build the command
    local cmd="${SCRIPT_DIR}/services/api/tests/run_tests.sh"
    
    # Add marker if specified
    if [[ -n "$backend_marker" ]]; then
        cmd="$cmd -m $backend_marker"
    fi
    
    # Add coverage if enabled
    if [[ "$coverage" = true ]]; then
        cmd="$cmd --html-report"
    else
        cmd="$cmd --no-coverage"
    fi
    
    # Make sure the script is executable
    chmod +x "$cmd"
    
    # Run the tests
    echo -e "${YELLOW}Running command:${NC} $cmd"
    $cmd
    return $?
}

# Function to run frontend tests
run_frontend_tests() {
    echo -e "\n${BOLD}${BLUE}=== Running Frontend Tests ===${NC}\n"
    
    # Check if frontend directory exists
    if [ ! -d "${SCRIPT_DIR}/frontend/testing-website" ]; then
        echo -e "${RED}Frontend directory not found at ${SCRIPT_DIR}/frontend/testing-website${NC}"
        return 1
    fi
    
    # Navigate to frontend directory
    cd "${SCRIPT_DIR}/frontend/testing-website"
    
    # Build the command
    local cmd="npm test"
    
    # Add coverage if enabled
    if [[ "$coverage" = true ]]; then
        cmd="$cmd -- --coverage --watchAll=false"
    else
        cmd="$cmd -- --watchAll=false"
    fi
    
    # Run the tests
    echo -e "${YELLOW}Running command:${NC} $cmd"
    eval $cmd
    return $?
}

# Track overall success
overall_success=true

# Run backend tests if requested
if [[ "$run_backend" = true ]]; then
    run_backend_tests
    if [[ $? -ne 0 ]]; then
        overall_success=false
        echo -e "${RED}Backend tests failed!${NC}"
    else
        echo -e "${GREEN}Backend tests passed!${NC}"
    fi
fi

# Run frontend tests if requested
if [[ "$run_frontend" = true ]]; then
    run_frontend_tests
    if [[ $? -ne 0 ]]; then
        overall_success=false
        echo -e "${RED}Frontend tests failed!${NC}"
    else
        echo -e "${GREEN}Frontend tests passed!${NC}"
    fi
fi

# Print summary
echo -e "\n${BOLD}${BLUE}=== Test Summary ===${NC}"
if [[ "$overall_success" = true ]]; then
    echo -e "${GREEN}All tests passed successfully!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed. Please check the output above for details.${NC}"
    exit 1
fi