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
        ("aripiprazole", 10),
        ("risperidone", 5),
        ("methylphenidate", 3)
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
    # The new format is a list of dicts with value and count
    med_values = [m["value"] for m in data["medication"]]
    assert "aripiprazole" in med_values
    assert "risperidone" in med_values
    assert "methylphenidate" in med_values
    
    # Verify that the execute method was called with the expected SQL
    assert mock_cursor.execute.call_count == 2
    
    # Check the first call (medications)
    call1 = mock_cursor.execute.call_args_list[0]
    sql1 = call1[0][0]
    assert 'SELECT LOWER(treatment_name)' in sql1
    assert "FROM jim_data.data_embedded" in sql1
    
    # Check the second call (symptoms)
    call2 = mock_cursor.execute.call_args_list[1]
    sql2 = call2[0][0]
    assert 'SELECT LOWER(primary_outcome_area)' in sql2
    assert "FROM jim_data.data_embedded" in sql2

@pytest.mark.api
def test_initial_results_endpoint(client, mock_db_connection):
    """Test the /api/initial-results endpoint returns properly formatted data."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response with all 38 columns from new schema
    mock_cursor.fetchall.return_value = [
        # title, pub_date, pmid, authors, url, Treatment name, Duration, 
        # Primary Outcome Area, Primary Outcome Measures, distance, abstract,
        # First Author, Date of publication, Study type, Sample Size, M:F Ratio,
        # Age range/mean, Medication/Treatment Dose Range, Results: Primary measure,
        # Secondary Outcome Area, Secondary Outcome Measures, Tolerability/Side Effects,
        # Safety, Drop Our Rate, Race/Ethnicity Percentages, Notes,
        # Sequence Generation, Allocation Concealment, Outcome Assessors Blinding,
        # Clinician and Participant Blinding, Incomplete outcome data, 
        # Selective outcome reporting, Notes on Biases, ai, age_min, age_max, 
        # males_in_study, females_in_study, Journal, Commercial Affiliation
        (
            "Test Study Title", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome Area", "Test Measures",
            0, "Test Abstract", "First Author", "2023-01-01", "RCT", "100", "1:1",
            "10-15", "10mg", "Improved", "Secondary Area", "Secondary Measures",
            "Mild", "Safe", "5%", "50% White", "Notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "No biases", "0.8",
            10, 15, 50, 50, "Test Journal", "Test Affiliation"
        ),
        (
            "Another Study", pub_date, 67890, "Another Author", "http://another.url",
            "Another Medication", "8 weeks", "Another Outcome", "Another Measures",
            0, "Another Abstract", "Second Author", "2023-02-01", "Open Label", "50",
            "2:1", "8-12", "15mg", "No change", "Another Secondary", "Other Measures",
            "Moderate", "Safe", "10%", "60% White", "Other notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "Some biases", "0.7",
            8, 12, 33, 17, "Another Journal", "Another Affiliation"
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
    assert "FROM jim_data.data_embedded" in sql_call
    assert "ORDER BY pub_date DESC" in sql_call

@pytest.mark.api
def test_search_endpoint_with_query(client, mock_db_connection, mock_sentence_model):
    """Test the /api/search endpoint with a query parameter."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response with all 38 columns from new schema
    mock_cursor.fetchall.return_value = [
        (
            "Test Study Title", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome Area", "Test Measures",
            0.2, "Test Abstract", "First Author", "2023-01-01", "RCT", "100", "1:1",
            "10-15", "10mg", "Improved", "Secondary Area", "Secondary Measures",
            "Mild", "Safe", "5%", "50% White", "Notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "No biases", "0.8",
            10, 15, 50, 50, "Test Journal", "Test Affiliation"
        ),
        (
            "Another Study", pub_date, 67890, "Another Author", "http://another.url",
            "Another Medication", "8 weeks", "Another Outcome", "Another Measures",
            0.3, "Another Abstract", "Second Author", "2023-02-01", "Open Label", "50",
            "2:1", "8-12", "15mg", "No change", "Another Secondary", "Other Measures",
            "Moderate", "Safe", "10%", "60% White", "Other notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "Some biases", "0.7",
            8, 12, 33, 17, "Another Journal", "Another Affiliation"
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
    assert "FROM jim_data.data_embedded" in sql_call
    assert "ORDER BY distance ASC" in sql_call

@pytest.mark.api
def test_search_endpoint_with_filters(client, mock_db_connection, mock_sentence_model):
    """Test the /api/search endpoint with filter parameters."""
    mock_conn, mock_cursor = mock_db_connection
    
    # Sample publication date for testing
    from datetime import date
    pub_date = date(2023, 1, 1)
    
    # Mock the database response with all 38 columns from new schema
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
    
    # Mock the database response for initial results with all 38 columns
    mock_cursor.fetchall.return_value = [
        (
            "Initial Result", pub_date, 12345, "Test Author", "http://test.url",
            "Test Medication", "6 weeks", "Test Outcome", "Test Measures",
            0, "Test Abstract", "First Author", "2023-01-01", "RCT", "100", "1:1",
            "10-15", "10mg", "Improved", "Secondary Area", "Secondary Measures",
            "Mild", "Safe", "5%", "50% White", "Notes", "Low risk", "Low risk",
            "Low risk", "Low risk", "Low risk", "Low risk", "No biases", "0.8",
            10, 15, 50, 50, "Test Journal", "Test Affiliation"
        )
    ]
    
    # Use a patch to mock the cache behavior
    with patch('app.search_cache.get') as mock_cache_get:
        # Set up the mock to return different values based on the key
        def side_effect(key):
            if key == 'initial_results_200':
                return None  # First call for initial_results returns None
            return None
        
        mock_cache_get.side_effect = side_effect
        
        # Call the API endpoint with an empty query
        response = client.get('/api/search')
        
        # Check that the response is successful
        assert response.status_code == 200
        
        # Verify it checked for initial_results in the cache
        mock_cache_get.assert_any_call('initial_results_200')