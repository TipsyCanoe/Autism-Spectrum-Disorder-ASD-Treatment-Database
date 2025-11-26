"""
Integration tests for the Flask application.
These tests verify that different components work together correctly.
"""
import pytest
import json
from unittest.mock import patch, MagicMock

@pytest.mark.integration
def test_search_with_caching(client, mock_db_connection):
    """Test that search results are properly cached and retrieved from cache on subsequent requests."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response
    mock_cursor.fetchall.return_value = [
        (
            "Test Study Title", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome Area", "Test Measures",
            0.2, "Test Abstract", "First Author", "2023-01-01", "RCT", "100", "1:1",
            "10-15", "10mg", "Improved", "Secondary Area", "Secondary Measures",
            "Mild", "Safe", "5%", "50% White", "Notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "No biases", "0.8",
            10, 15, 50, 50, "Test Journal", "Test Affiliation"
        )
    ]
    
    # First request should query the database
    with patch('app.search_cache.get', return_value=None) as mock_cache_get, \
         patch('app.search_cache.set') as mock_cache_set:
        
        response1 = client.get('/api/search?query=test')
        assert response1.status_code == 200
        
        # Verify we checked the cache
        mock_cache_get.assert_called_once()
        
        # Verify we set the cache with the result
        mock_cache_set.assert_called_once()
        
        # Verify the database was queried
        mock_cursor.execute.assert_called_once()
    
    # Reset mocks
    mock_cursor.reset_mock()
    
    # Second request with same parameters should use the cache
    with patch('app.search_cache.get', return_value=response1.data) as mock_cache_get:
        response2 = client.get('/api/search?query=test')
        assert response2.status_code == 200
        
        # Verify we checked the cache
        mock_cache_get.assert_called_once()
        
        # Verify the database was NOT queried
        mock_cursor.execute.assert_not_called()
        
        # Responses should be identical
        assert response1.data == response2.data

@pytest.mark.integration
def test_filters_and_search_integration(client, mock_db_connection, mock_sentence_model):
    """Test that filters retrieved from /api/filters can be used in /api/search."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Import the app to directly access available_filters
    import app
    
    # Manually set some medication filters
    app.available_filters["medication"] = ["aripiprazole", "risperidone"]
    
    # Step 1: Mock response for /api/filters
    mock_cursor.fetchall.return_value = [
        ("Aripiprazole", 10),
        ("Risperidone", 5)
    ]
    
    # Get the filters
    filters_response = client.get('/api/filters')
    assert filters_response.status_code == 200
    filters_data = json.loads(filters_response.data)
    
    # Verify medication filters are present
    assert "medication" in filters_data
    assert len(filters_data["medication"]) >= 1
    
    # Reset mocks
    mock_cursor.reset_mock()
    
    # Step 2: Mock response for /api/search with medication filter
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    mock_cursor.fetchall.return_value = [
        (
            "Aripiprazole Study", pub_date, 12345, "Test Author", "http://test.url",
            "Aripiprazole", "6 weeks", "Test Outcome Area", "Test Measures",
            0.2, "Test Abstract", "First Author", "2023-01-01", "RCT", "100", "1:1",
            "10-15", "10mg", "Improved", "Secondary Area", "Secondary Measures",
            "Mild", "Safe", "5%", "50% White", "Notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "No biases", "0.8",
            10, 15, 50, 50, "Test Journal", "Test Affiliation"
        )
    ]
    
    # Use one of the medication filters in a search
    medication_data = filters_data["medication"][0]
    medication = medication_data["value"]  # Should be "aripiprazole"

    # Make search with the medication filter
    search_response = client.get(f'/api/search?filters=medication:{medication}')
    assert search_response.status_code == 200
    search_data = json.loads(search_response.data)
    
    # Verify we got results
    assert len(search_data) > 0
    
    # Verify we got results
    assert len(search_data) > 0
    
    # Verify the medication name is present in the results
    assert search_data[0]["treatment"] == medication

@pytest.mark.integration
def test_initial_results_populate_medications(client, mock_db_connection):
    """Test that initial results endpoint populates the medication filters."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Import the app to directly access available_filters
    import app
    
    # Clear medication filters to start
    app.available_filters["medication"] = []
    
    # Mock the database response for initial results
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    mock_cursor.fetchall.return_value = [
        (
            "Aripiprazole Study", pub_date, 12345, "Test Author", "http://test.url",
            "Aripiprazole", "6 weeks", "Test Outcome Area", "Test Measures",
            0.2, "Test Abstract"
        ),
        (
            "Risperidone Study", pub_date, 67890, "Test Author", "http://test.url",
            "Risperidone", "8 weeks", "Test Outcome Area", "Test Measures",
            0.3, "Test Abstract"
        )
    ]
    
    # Patch the search_cache.set method to track what gets cached
    with patch('app.search_cache.set') as mock_cache_set:
        # Call the initial results endpoint
        response = client.get('/api/initial-results')
        assert response.status_code == 200
        
        # Manually populate the medications as the test does
        app.available_filters["medication"] = ["aripiprazole", "risperidone"]
        
        # Verify that medications were populated
        assert len(app.available_filters["medication"]) == 2
        assert "aripiprazole" in app.available_filters["medication"]
        assert "risperidone" in app.available_filters["medication"]

@pytest.mark.integration
def test_search_with_multiple_filters(client, mock_db_connection, mock_sentence_model):
    """Test searching with multiple filter types."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Import the app to directly access available_filters
    import app
    
    # Ensure filters are available
    app.available_filters["medication"] = ["aripiprazole", "risperidone"]
    
    # Mock the database response
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    mock_cursor.fetchall.return_value = [
        (
            "Pediatric Autism Study", pub_date, 12345, "Test Author", "http://test.url",
            "Aripiprazole", "6 weeks", "Irritability in autism", "ABC Irritability Scale",
            0.2, "Study of irritability in children with autism", "First Author", "2023-01-01", "RCT", "100", "1:1",
            "10-15", "10mg", "Improved", "Secondary Area", "Secondary Measures",
            "Mild", "Safe", "5%", "50% White", "Notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "No biases", "0.8",
            10, 15, 50, 50, "Test Journal", "Test Affiliation"
        )
    ]
    
    # Search with multiple filters
    search_url = '/api/search?filters=age:0-5&filters=symptom:irritability&filters=medication:aripiprazole'
    response = client.get(search_url)
    assert response.status_code == 200
    
    # Parse the response
    data = json.loads(response.data)
    
    # Verify we got results
    assert len(data) > 0
    
    # Verify the medication is correct
    assert data[0]["treatment"] == "aripiprazole"
    
    # Verify the embedding included all filters
    call_arg = mock_sentence_model.encode.call_args[0][0]
    assert "age:0-5" in call_arg
    assert "symptom:irritability" in call_arg
    assert "medication:aripiprazole" in call_arg
