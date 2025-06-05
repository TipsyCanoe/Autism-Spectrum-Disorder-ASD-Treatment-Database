import pytest
import sys
import os
import json
from unittest.mock import patch, MagicMock

# Add the parent directory to path so we can import app
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app import app as flask_app

@pytest.fixture
def app():
    """Create and configure a Flask app for testing."""
    # Configure app for testing
    flask_app.config.update({
        "TESTING": True,
    })
    return flask_app

@pytest.fixture
def client(app):
    """A test client for the app."""
    return app.test_client()

@pytest.fixture
def mock_db_connection():
    """Mock database connection for tests."""
    mock_conn = MagicMock()
    mock_cursor = MagicMock()
    mock_conn.cursor.return_value = mock_cursor
    
    # Mock the execute and fetchall methods
    mock_cursor.execute = MagicMock()
    mock_cursor.fetchall = MagicMock()
    
    with patch('app.get_db_connection', return_value=mock_conn):
        yield mock_conn, mock_cursor

@pytest.fixture
def sample_study():
    """Return a sample study dictionary for testing filter matching functions."""
    return {
        "Study Title": "Effect of aripiprazole on autism symptoms in children",
        "Primary Outcome Area": "irritability and hyperactivity",
        "Secondary Outcome Area": "social functioning",
        "Results: Primary measure": "Reduced irritability scores",
        "Results: Secondary Measures": "Improved social interaction",
        "Tolerability/Side Effects": "Mild sedation",
        "Safety": "Generally well tolerated",
        "Age Range": "8-17 / 12.5",
        "M:F": "3:1"
    }

@pytest.fixture
def mock_sentence_model():
    """Mock the sentence transformer model."""
    import numpy as np
    mock_model = MagicMock()
    mock_model.encode.return_value = np.array([0.1, 0.2, 0.3, 0.4])  # Return numpy array instead of list
    
    with patch('app.sentence_model', mock_model):
        yield mock_model
