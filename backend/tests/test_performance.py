"""
Performance tests for the Flask API.
These tests measure response times and validate that caching improves performance.
"""
import pytest
import time
import json
from flask import jsonify
from unittest.mock import patch, MagicMock

@pytest.mark.performance
def test_search_endpoint_performance(client, mock_db_connection, mock_sentence_model):
    """Test the performance improvement from caching for the search endpoint."""
    mock_conn, mock_cursor = mock_db_connection

    # Mock the database response
    from datetime import date
    pub_date = date(2023, 1, 1)

    mock_cursor.fetchall.return_value = [
        (
            "Test Study", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome", "Test Measures",
            0.2, "Test Abstract"
        )
    ]

    # Import the app to directly access and reset caches
    import app
    app.search_cache = app.SimpleCache(max_size=10)

    # Use a specific query to ensure cache key consistency
    query = "autism children"

    # First request (no cache)
    with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db_first:
        start_time = time.time()
        response1 = client.get(f'/api/search?query={query}')
        first_request_time = time.time() - start_time
        assert response1.status_code == 200
        
        # Verify database was accessed
        assert mock_get_db_first.called

    # Second request (with cache)
    with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db_second:
        start_time = time.time()
        response2 = client.get(f'/api/search?query={query}')
        second_request_time = time.time() - start_time
        assert response2.status_code == 200
        
        # Verify database was NOT accessed (used cache)
        assert not mock_get_db_second.called

    # Log the performance improvement
    improvement_factor = first_request_time / second_request_time if second_request_time > 0 else float('inf')
    print(f"\nSearch endpoint performance:")
    print(f"  First request: {first_request_time:.6f} seconds")
    print(f"  Cached request: {second_request_time:.6f} seconds")
    print(f"  Improvement factor: {improvement_factor:.2f}x faster with caching")

@pytest.mark.performance
def test_initial_results_performance(client, mock_db_connection):
    """Test the performance of the initial results endpoint with caching."""
    mock_conn, mock_cursor = mock_db_connection

    # Mock the database response
    from datetime import date
    pub_date = date(2023, 1, 1)

    mock_cursor.fetchall.return_value = [
        (
            "Test Study", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome", "Test Measures",
            0, "Test Abstract"
        )
    ]

    # Import the app to directly access and reset caches
    import app
    app.search_cache = app.SimpleCache(max_size=10)

    # First request (no cache)
    with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db_first:
        start_time = time.time()
        response1 = client.get('/api/initial-results')
        first_request_time = time.time() - start_time
        assert response1.status_code == 200
        
        # Verify database was accessed
        assert mock_get_db_first.called

    # Second request (with cache)
    with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db_second:
        start_time = time.time()
        response2 = client.get('/api/initial-results')
        second_request_time = time.time() - start_time
        assert response2.status_code == 200
        
        # Verify database was NOT accessed (used cache)
        assert not mock_get_db_second.called

    # Log the performance improvement
    improvement_factor = first_request_time / second_request_time if second_request_time > 0 else float('inf')
    print(f"\nInitial results endpoint performance:")
    print(f"  First request: {first_request_time:.6f} seconds")
    print(f"  Cached request: {second_request_time:.6f} seconds")
    print(f"  Improvement factor: {improvement_factor:.2f}x faster with caching")

@pytest.mark.performance
def test_load_simulation(client, mock_db_connection, mock_sentence_model):
    """
    Simulate multiple users making requests to test cache efficiency
    under load.
    """
    mock_conn, mock_cursor = mock_db_connection

    # Mock the database response
    from datetime import date
    pub_date = date(2023, 1, 1)

    mock_cursor.fetchall.return_value = [
        (
            "Test Study", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome", "Test Measures",
            0.2, "Test Abstract"
        )
    ]

    # Import the app to directly access and reset caches
    import app
    app.search_cache = app.SimpleCache(max_size=10)
    app.filter_cache = app.SimpleCache(max_size=10)

    # Define common search queries
    queries = [
        "autism",
        "adhd",
        "hyperactivity",
        "irritability",
        "aripiprazole"
    ]

    # Track database access counts
    db_access_count = 0

    # Simulate 10 requests (simulating multiple users)
    for i in range(10):
        # Reset mocks for clean tracking
        mock_conn.reset_mock()
        
        # Alternate between different endpoints
        if i % 3 == 0:
            # Get filters
            with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db:
                response = client.get('/api/filters')
                if mock_get_db.called:
                    db_access_count += 1
                
        elif i % 3 == 1:
            # Get initial results
            with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db:
                response = client.get('/api/initial-results')
                if mock_get_db.called:
                    db_access_count += 1
                
        else:
            # Search with different queries
            query = queries[i % len(queries)]
            with patch('app.get_db_connection', return_value=mock_conn) as mock_get_db:
                response = client.get(f'/api/search?query={query}')
                if mock_get_db.called:
                    db_access_count += 1

        assert response.status_code == 200

    # For 10 requests, we should see fewer than 10 database accesses due to caching
    print(f"\nLoad simulation results:")
    print(f"  Total requests: 10")
    print(f"  Database accesses: {db_access_count}")
    print(f"  Cache efficiency: {(10 - db_access_count) / 10 * 100:.1f}%")
    
    # If caching is working, we should have fewer DB accesses than requests
    assert db_access_count < 10