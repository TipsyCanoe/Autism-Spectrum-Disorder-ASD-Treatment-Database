import React from 'react';

// Mock the App component (top level)
jest.mock('./app/App', () => () => <div data-testid="app-mock">AppMock</div>);

// Define mocks that will be used by the react-dom/client mock (top level)
const mockRender = jest.fn();
const mockCreateRoot = jest.fn(() => ({ // Default implementation
  render: mockRender,
}));

// Mock react-dom/client (top level)
// This tells Jest to replace the actual 'react-dom/client' module with this mock object.
jest.mock('react-dom/client', () => ({
  // Spread properties from the actual module first
  ...jest.requireActual('react-dom/client'),
  // Override createRoot with our top-level mock function
  createRoot: mockCreateRoot,
}));

describe('src/index.js', () => {
  let rootElement;

  beforeEach(() => {
    // Reset modules to ensure index.js runs fresh for each test,
    // as it has side effects (calling ReactDOM.createRoot().render()).
    jest.resetModules();

    // Clear call history and reset implementations for mocks before each test
    mockRender.mockClear();
    mockCreateRoot.mockClear();
    // Ensure mockCreateRoot returns an object with the current mockRender for each test run
    mockCreateRoot.mockImplementation(() => ({
      render: mockRender,
    }));

    // Create a mock root element and append it to the body
    rootElement = document.createElement('div');
    rootElement.id = 'root';
    document.body.appendChild(rootElement);

    // Spy on document.getElementById to control its behavior for this test
    jest.spyOn(document, 'getElementById').mockImplementation((id) => {
      if (id === 'root') {
        return rootElement;
      }
      return null;
    });
  });

  afterEach(() => {
    // Clean up the root element from the document body
    if (rootElement && rootElement.parentNode) {
      rootElement.parentNode.removeChild(rootElement);
    }
    // Restore all mocks and spies (e.g., document.getElementById)
    jest.restoreAllMocks();
  });

  test('should call ReactDOM.createRoot with the root element and render the App component in StrictMode', () => {
    // Dynamically import/require index.js to execute its code now that mocks are set up.
    // jest.resetModules() in beforeEach ensures this is a fresh execution.
    require('../src/index.js');

    // Check if document.getElementById was called with 'root'
    expect(document.getElementById).toHaveBeenCalledWith('root');

    // Check if ReactDOM.createRoot (our mock) was called with the root element
    expect(mockCreateRoot).toHaveBeenCalledWith(rootElement);

    // Check if the render method (our mockRender) was called once
    expect(mockRender).toHaveBeenCalledTimes(1);

    // Check if render was called with <React.StrictMode><App /></React.StrictMode>
    const renderCallArgs = mockRender.mock.calls[0][0];
    expect(renderCallArgs.type).toBe(React.StrictMode);

    // Verify that the child of StrictMode is the mocked App component
    const AppMockComponent = jest.requireMock('./app/App');
    expect(renderCallArgs.props.children.type).toBe(AppMockComponent);
  });
});
