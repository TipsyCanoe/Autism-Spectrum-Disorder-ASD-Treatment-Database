import React from "react";

const FilterSection = ({
  title,
  category,
  options,
  selectedOptions,
  onFilterChange,
  isExpanded,
  onToggle,
}) => {
  return (
    <div className="mb-6 border-t border-gray-200 pt-4">
      <button
        className="flex justify-between items-center w-full text-left font-medium text-gray-700 mb-2"
        onClick={onToggle}
        aria-expanded={isExpanded} // Add aria-expanded for accessibility
      >
        <span>{title}</span>
        <span>{isExpanded ? "−" : "+"}</span>
      </button>

      {isExpanded && (
        <div className="ml-2 space-y-2">
          {options.map((option) => (
            <label key={option.value} className="flex items-center">
              <input
                type="checkbox"
                className="mr-2"
                checked={selectedOptions.includes(
                  `${category}:${option.value}`
                )}
                onChange={() => onFilterChange(category, option.value)}
              />
              <span>{option.label}</span>
            </label>
          ))}
        </div>
      )}
    </div>
  );
};

export default FilterSection;
