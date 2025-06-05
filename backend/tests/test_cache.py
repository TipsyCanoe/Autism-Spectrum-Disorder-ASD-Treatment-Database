"""
Tests for the SimpleCache implementation used in the Flask application.
"""
import pytest
from app import SimpleCache

@pytest.mark.unit
def test_simple_cache_get_set():
    """Test basic get and set functionality of SimpleCache."""
    cache = SimpleCache(max_size=3)
    
    # Test setting and getting a value
    cache.set("key1", "value1")
    assert cache.get("key1") == "value1"
    
    # Test getting a non-existent key
    assert cache.get("nonexistent") is None
    
    # Test overwriting a value
    cache.set("key1", "updated_value")
    assert cache.get("key1") == "updated_value"

@pytest.mark.unit
def test_simple_cache_max_size():
    """Test that SimpleCache respects max_size and removes oldest items."""
    cache = SimpleCache(max_size=2)
    
    # Add two items to fill the cache
    cache.set("key1", "value1")
    cache.set("key2", "value2")
    
    # Verify both items are in the cache
    assert cache.get("key1") == "value1"
    assert cache.get("key2") == "value2"
    
    # Add a third item, which should evict the first item
    cache.set("key3", "value3")
    
    # Verify the first item is gone and the other two remain
    assert cache.get("key1") is None
    assert cache.get("key2") == "value2"
    assert cache.get("key3") == "value3"

@pytest.mark.unit
def test_simple_cache_with_complex_values():
    """Test that SimpleCache can handle complex values like dictionaries and lists."""
    cache = SimpleCache()
    
    # Test with a dictionary
    dict_value = {"name": "test", "values": [1, 2, 3]}
    cache.set("dict_key", dict_value)
    assert cache.get("dict_key") == dict_value
    
    # Test with a list
    list_value = [1, "two", {"three": 3}]
    cache.set("list_key", list_value)
    assert cache.get("list_key") == list_value
    
    # Test with None
    cache.set("none_key", None)
    assert cache.get("none_key") is None
    # Make sure we can distinguish between a key with None and a missing key
    assert "none_key" in cache.cache
    assert "missing_key" not in cache.cache

@pytest.mark.unit
def test_cache_size_zero():
    """Test SimpleCache behavior with max_size of 0 (edge case)."""
    # SimpleCache should enforce a minimum size of 1
    cache = SimpleCache(max_size=0)
    
    # Should still allow at least one item because the implementation
    # enforces a minimum size of 1
    cache.set("key1", "value1")
    assert cache.get("key1") == "value1"
    
    # Adding a second item should replace the first due to max_size=1
    cache.set("key2", "value2")
    
    # Verify only the second item remains
    assert cache.get("key1") is None
    assert cache.get("key2") == "value2"
