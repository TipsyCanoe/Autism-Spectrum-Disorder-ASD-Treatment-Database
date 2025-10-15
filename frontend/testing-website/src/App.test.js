// filepath: /home/coleoliva/senior-proj/frontend/testing-website/src/App.test.js
import { render, screen, within } from '@testing-library/react'; // Import within
import App from './app/App'; // Assuming App.jsx is in src/app/App.jsx

describe('App Component Rendering and Initial Route', () => {
  test('renders the App without crashing and displays the NavBar', () => {
    render(<App />);

    // Verify NavBar title is rendered (already present)
    const navTitle = screen.getByText(/STAR/i);
    expect(navTitle).toBeInTheDocument();

    // Find the navigation bar element first
    const navBar = screen.getByRole('navigation'); // Assuming your NavBar has a role="navigation"

    // Verify NavBar links are present *within* the navBar
    expect(within(navBar).getByRole('link', { name: /Home/i })).toBeInTheDocument();
    expect(within(navBar).getByRole('link', { name: /Search/i })).toBeInTheDocument(); // Use within here
    expect(within(navBar).getByRole('link', { name: /FAQ/i })).toBeInTheDocument();
    expect(within(navBar).getByRole('link', { name: /About/i })).toBeInTheDocument();
  });

  test('renders the HomePage component on the default route', () => {
    render(<App />);

    // Verify a key element from HomePage is rendered
    const welcomeHeading = screen.getByRole('heading', {
      name: /Welcome to the Autism Resources Database/i,
    });
    expect(welcomeHeading).toBeInTheDocument();

    // This test can check for the button on the HomePage.
    // To be more specific, you could find the main content area first.
    expect(screen.getByRole('button', { name: /Search/i })).toBeInTheDocument();
  });
});
