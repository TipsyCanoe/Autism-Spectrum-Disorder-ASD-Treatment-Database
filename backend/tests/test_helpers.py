import pytest
import sys
import os
import json
from unittest.mock import patch, MagicMock

# Add the parent directory to path so we can import app functions
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app import (
    extract_age_from_range,
    matches_age_filter,
    matches_gender_filter,
    matches_symptom_filter,
    matches_text_query,
    SimpleCache
)

@pytest.mark.unit
def test_extract_age_from_range():
    """Test the extraction of age ranges from string."""
    # Test valid ranges
    assert extract_age_from_range("8-17 / 12.5") == (8, 17)
    assert extract_age_from_range("6-12 / 9") == (6, 12)
    assert extract_age_from_range("18-65 / 34.2") == (18, 65)
    
    # Test edge cases
    assert extract_age_from_range("") is None
    assert extract_age_from_range(None) is None
    assert extract_age_from_range("Not a range") is None
    assert extract_age_from_range("5 years") is None

@pytest.mark.unit
def test_matches_age_filter():
    """Test matching studies to age filters."""
    # Study ages that should match specific ranges
    assert matches_age_filter("3-5 / 4", ["0-5"]) is True
    assert matches_age_filter("10-15 / 12", ["6-12"]) is True
    assert matches_age_filter("15-16 / 15.5", ["13-17"]) is True
    assert matches_age_filter("20-30 / 25", ["18-25"]) is True  # Overlaps
    assert matches_age_filter("50-60 / 55", ["26-64"]) is True
    assert matches_age_filter("70-80 / 75", ["65+"]) is True
    
    # Study ages that should not match specific ranges
    assert matches_age_filter("6-10 / 8", ["0-5"]) is False
    assert matches_age_filter("18-25 / 21", ["13-17"]) is False
    
    # Invalid age ranges
    assert matches_age_filter("", ["0-5"]) is False
    assert matches_age_filter(None, ["0-5"]) is False
    assert matches_age_filter("Not a range", ["0-5"]) is False

@pytest.mark.unit
def test_matches_gender_filter():
    """Test matching studies to gender filters."""
    # Valid gender ratios
    assert matches_gender_filter("3:1", ["male"]) is True
    assert matches_gender_filter("1:3", ["female"]) is True
    assert matches_gender_filter("1:1", ["male", "female"]) is True
    
    # Invalid inputs
    assert matches_gender_filter("", ["male"]) is False
    assert matches_gender_filter(None, ["male"]) is False
    
    # TODO: Add tests for nonbinary if implementation supports it

@pytest.mark.unit
def test_matches_symptom_filter(sample_study):
    """Test matching studies to symptom filters."""
    # The sample study has "irritability and hyperactivity" in Primary Outcome Area
    assert matches_symptom_filter(sample_study, ["irritability"]) is True
    assert matches_symptom_filter(sample_study, ["hyperactivity"]) is True
    
    # The sample study has "social functioning" in Secondary Outcome Area
    assert matches_symptom_filter(sample_study, ["social"]) is True
    
    # The sample study doesn't mention anxiety
    assert matches_symptom_filter(sample_study, ["anxiety-reactivity"]) is False
    
    # Test with empty study
    assert matches_symptom_filter({}, ["irritability"]) is False

@pytest.mark.unit
def test_matches_text_query(sample_study):
    """Test matching studies to text queries."""
    # Simple matches
    assert matches_text_query(sample_study, "autism") is True
    assert matches_text_query(sample_study, "children") is True
    assert matches_text_query(sample_study, "aripiprazole") is True
    
    # Multi-word matches (all terms must be present)
    assert matches_text_query(sample_study, "autism children") is True
    assert matches_text_query(sample_study, "aripiprazole autism") is True
    
    # Non-matches
    assert matches_text_query(sample_study, "risperidone") is False
    assert matches_text_query(sample_study, "adults depression") is False
    
    # Empty query always matches
    assert matches_text_query(sample_study, "") is True
    
    # Empty study
    assert matches_text_query({}, "autism") is False

@pytest.mark.unit
def test_simple_cache():
    """Test the SimpleCache implementation."""
    # Create a cache with max size 2
    cache = SimpleCache(max_size=2)
    
    # Test basic get/set
    assert cache.get("key1") is None  # Key doesn't exist yet
    
    cache.set("key1", "value1")
    assert cache.get("key1") == "value1"
    
    # Test eviction when full
    cache.set("key2", "value2")
    cache.set("key3", "value3")  # This should evict key1
    
    assert cache.get("key1") is None  # key1 should be evicted
    assert cache.get("key2") == "value2"
    assert cache.get("key3") == "value3"
    
    # Test overwriting existing key
    cache.set("key2", "new_value2")
    assert cache.get("key2") == "new_value2"
