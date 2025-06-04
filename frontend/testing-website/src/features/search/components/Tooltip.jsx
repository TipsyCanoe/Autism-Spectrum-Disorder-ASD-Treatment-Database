import PropTypes from 'prop-types';
import { useState } from 'react';

const Tooltip = ({ text, children }) => {
  const [showTooltip, setShowTooltip] = useState(false);

  return (
    <div 
      className="relative inline-block"
      onMouseEnter={() => setShowTooltip(true)}
      onMouseLeave={() => setShowTooltip(false)}
    >
      {children}
      {showTooltip && text && (
        <div 
          className="absolute z-10 px-3 py-2 text-sm font-medium text-white bg-gray-900 rounded-lg shadow-sm tooltip bottom-full left-1/2 transform -translate-x-1/2 mb-2 whitespace-normal max-w-xs"
          style={{ minWidth: '100px' }} // Ensure tooltip has some width
        >
          {text}
          <div className="tooltip-arrow" data-popper-arrow></div>
        </div>
      )}
    </div>
  );
};

Tooltip.propTypes = {
  text: PropTypes.string,
  children: PropTypes.node.isRequired,
};

export default Tooltip;
