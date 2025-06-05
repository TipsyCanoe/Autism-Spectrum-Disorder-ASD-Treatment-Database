"""
Tests for API endpoints in the Flask application.
This file contains tests that make actual HTTP requests to the API endpoints.
"""
import pytest
import json
from unittest.mock import patch, MagicMock

@pytest.mark.api
def test_get_filters_endpoint(client, mock_db_connection):
    """Test the /api/filters endpoint returns the expected filter structure."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Mock the database response for medications
    mock_cursor.fetchall.return_value = [
        ("Aripiprazole",),
        ("Risperidone",),
        ("Methylphenidate",)
    ]
    
    # Call the API endpoint
    response = client.get('/api/filters')
    
    # Check that the response is successful
    assert response.status_code == 200
    
    # Parse the response data
    data = json.loads(response.data)
    
    # Check that the response contains the expected filter categories
    assert "age" in data
    assert "symptom" in data
    assert "gender" in data
    assert "medication" in data
    
    # Check that the medication list contains the mocked values
    assert "aripiprazole" in data["medication"]
    assert "risperidone" in data["medication"]
    assert "methylphenidate" in data["medication"]
    
    # Verify that the execute method was called with the expected SQL
    mock_cursor.execute.assert_called_once()
    sql_call = mock_cursor.execute.call_args[0][0]
    assert "SELECT DISTINCT treatment_name" in sql_call
    assert "FROM updated_treatment_data.treatments" in sql_call

@pytest.mark.api
def test_initial_results_endpoint(client, mock_db_connection):
    """Test the /api/initial-results endpoint returns properly formatted data."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response
    mock_cursor.fetchall.return_value = [
        # title, pub_date, pmid, authors, url, treatment_name, duration, 
        # primary_outcome_area, primary_outcome_measures, distance, abstract
        (
            "Test Study Title", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome Area", "Test Measures",
            0, "Test Abstract"
        ),
        (
            "Another Study", pub_date, 67890, "Another Author", "http://another.url",
            "Another Medication", "8 weeks", "Another Outcome", "Another Measures",
            0, "Another Abstract"
        )
    ]
    
    # Call the API endpoint
    response = client.get('/api/initial-results')
    
    # Check that the response is successful
    assert response.status_code == 200
    
    # Parse the response data
    data = json.loads(response.data)
    
    # Check the structure of the response
    assert isinstance(data, list)
    assert len(data) > 0
    
    # Check the first entry in the response
    first_entry = data[0]
    assert "treatment" in first_entry
    assert "studies" in first_entry
    assert isinstance(first_entry["studies"], list)
    
    # Check a study entry
    if first_entry["studies"]:
        first_study = first_entry["studies"][0]
        assert "Study Title" in first_study
        assert "PMID" in first_study
        assert "Publication Date" in first_study
        assert "Abstract" in first_study
    
    # Verify that the execute method was called with the expected SQL
    mock_cursor.execute.assert_called_once()
    sql_call = mock_cursor.execute.call_args[0][0]
    assert "FROM updated_treatment_data.semantic_paper_search_view" in sql_call
    assert "ORDER BY spv.pub_date DESC" in sql_call

@pytest.mark.api
def test_search_endpoint_with_query(client, mock_db_connection, mock_sentence_model):
    """Test the /api/search endpoint with a query parameter."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response
    mock_cursor.fetchall.return_value = [
        # Similar structure as initial_results test but with a distance value
        (
            "Test Study Title", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome Area", "Test Measures",
            0.2, "Test Abstract"
        ),
        (
            "Another Study", pub_date, 67890, "Another Author", "http://another.url",
            "Another Medication", "8 weeks", "Another Outcome", "Another Measures",
            0.3, "Another Abstract"
        )
    ]
    
    # Call the API endpoint with a query
    response = client.get('/api/search?query=autism')
    
    # Check that the response is successful
    assert response.status_code == 200
    
    # Parse the response data
    data = json.loads(response.data)
    
    # Check the structure of the response
    assert isinstance(data, list)
    assert len(data) > 0
    
    # Check the first entry in the response
    first_entry = data[0]
    assert "treatment" in first_entry
    assert "studies" in first_entry
    assert isinstance(first_entry["studies"], list)
    
    # Check a study entry
    if first_entry["studies"]:
        first_study = first_entry["studies"][0]
        assert "Study Title" in first_study
        assert "PMID" in first_study
        assert "Similarity Score" in first_study
        assert "Distance" in first_study
    
    # Verify the model was used to create embeddings
    mock_sentence_model.encode.assert_called_once()
    
    # Verify that the execute method was called with the expected SQL
    mock_cursor.execute.assert_called_once()
    sql_call = mock_cursor.execute.call_args[0][0]
    assert "FROM updated_treatment_data.semantic_paper_search_view" in sql_call
    assert "ORDER BY distance ASC" in sql_call

@pytest.mark.api
def test_search_endpoint_with_filters(client, mock_db_connection, mock_sentence_model):
    """Test the /api/search endpoint with filter parameters."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response
    mock_cursor.fetchall.return_value = [
        # Similar structure as search_with_query test
        (
            "Test Study Title", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome Area", "Test Measures",
            0.2, "Test Abstract"
        )
    ]
    
    # Call the API endpoint with filters
    response = client.get('/api/search?filters=age:0-5&filters=symptom:irritability')
    
    # Check that the response is successful
    assert response.status_code == 200
    
    # Parse the response data
    data = json.loads(response.data)
    
    # Verify the model was used to create embeddings with filter info
    mock_sentence_model.encode.assert_called_once()
    encode_arg = mock_sentence_model.encode.call_args[0][0]
    assert "age:0-5, symptom:irritability" in encode_arg
    
    # Verify that the execute method was called with the expected SQL
    mock_cursor.execute.assert_called_once()

@pytest.mark.api
def test_empty_search_uses_initial_results(client, mock_db_connection):
    """Test that an empty search redirects to initial results."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response for initial results
    mock_cursor.fetchall.return_value = [
        (
            "Initial Result", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome", "Test Measures",
            0, "Test Abstract"
        )
    ]
    
    # Use a patch to mock the cache behavior
    with patch('app.search_cache.get') as mock_cache_get:
        # Set up the mock to return different values based on the key
        def side_effect(key):
            if key == 'initial_results':
                return None  # First call for initial_results returns None
            return None
        
        mock_cache_get.side_effect = side_effect
        
        # Call the API endpoint with an empty query
        response = client.get('/api/search')
        
        # Check that the response is successful
        assert response.status_code == 200
        
        # Verify it checked for initial_results in the cache
        mock_cache_get.assert_any_call('initial_results')
