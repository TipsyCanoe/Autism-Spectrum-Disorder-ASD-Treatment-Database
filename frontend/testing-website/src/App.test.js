import { render, screen } from '@testing-library/react';
import App from './app/App';
import axios from 'axios';

// Mock axios so its import in searchService doesn't cause errors.
jest.mock('axios');

test('renders the App without crashing', () => {
  render(<App />);
  // Verify that a key element from your root layout (e.g., NavBar title "ARD") is rendered
  const navTitle = screen.getByText(/ARD/i);
  expect(navTitle).toBeInTheDocument();
});