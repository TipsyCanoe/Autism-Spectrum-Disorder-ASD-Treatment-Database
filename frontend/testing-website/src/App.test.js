// filepath: /home/coleoliva/senior-proj/frontend/testing-website/src/App.test.js
import { render, screen, within } from '@testing-library/react'; // Import within
import App from './app/App'; // Assuming App.jsx is in src/app/App.jsx

describe('App Component Rendering and Initial Route', () => {
  test('renders the App without crashing and displays the NavBar', () => {
    render(<App />);

    // Verify NavBar title is rendered - be more specific to avoid matching homepage content
    const navBar = screen.getByRole('navigation');
    const navTitle = within(navBar).getByText(/STAR/i);
    expect(navTitle).toBeInTheDocument();

    // Verify NavBar links are present *within* the navBar
    expect(within(navBar).getByRole('link', { name: /Home/i })).toBeInTheDocument();
    expect(within(navBar).getByRole('link', { name: /Search/i })).toBeInTheDocument();
    expect(within(navBar).getByRole('link', { name: /FAQ/i })).toBeInTheDocument();
    expect(within(navBar).getByRole('link', { name: /About/i })).toBeInTheDocument();
  });

  test('renders the HomePage component on the default route', () => {
    render(<App />);

    // Verify a key element from HomePage is rendered - updated to match new heading
    const welcomeHeading = screen.getByRole('heading', {
      name: /Welcome to STAR/i,
    });
    expect(welcomeHeading).toBeInTheDocument();

    // This test can check for the button on the HomePage.
    // To be more specific, you could find the main content area first.
    expect(screen.getByRole('button', { name: /Search/i })).toBeInTheDocument();
  });
});
