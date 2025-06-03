import '@testing-library/jest-dom';
import { fireEvent, render, screen } from '@testing-library/react';
import FilterSection from './FilterSection';

describe('FilterSection Component', () => {
  const mockOnFilterChange = jest.fn();
  const mockOnToggle = jest.fn();
  const defaultProps = {
    title: 'Test Category',
    category: 'testCategory',
    options: [
      { value: 'opt1', label: 'Option 1' },
      { value: 'opt2', label: 'Option 2' },
    ],
    selectedOptions: ['testCategory:opt1'],
    onFilterChange: mockOnFilterChange,
    isExpanded: true,
    onToggle: mockOnToggle,
  };

  beforeEach(() => {
    jest.clearAllMocks();
  });

  test('renders title correctly', () => {
    render(<FilterSection {...defaultProps} />);
    expect(screen.getByText('Test Category')).toBeInTheDocument();
  });

  test('renders expand/collapse button with correct aria-expanded state and text', () => {
    const { rerender } = render(<FilterSection {...defaultProps} isExpanded={true} />);
    const toggleButton = screen.getByRole('button', { name: /Test Category/i });
    expect(toggleButton).toHaveAttribute('aria-expanded', 'true');
    expect(screen.getByText('−')).toBeInTheDocument(); // Expanded indicator

    rerender(<FilterSection {...defaultProps} isExpanded={false} />);
    expect(toggleButton).toHaveAttribute('aria-expanded', 'false');
    expect(screen.getByText('+')).toBeInTheDocument(); // Collapsed indicator
  });

  test('calls onToggle when the toggle button is clicked', () => {
    render(<FilterSection {...defaultProps} />);
    const toggleButton = screen.getByRole('button', { name: /Test Category/i });
    fireEvent.click(toggleButton);
    expect(mockOnToggle).toHaveBeenCalledTimes(1);
  });

  test('renders options when isExpanded is true', () => {
    render(<FilterSection {...defaultProps} isExpanded={true} />);
    expect(screen.getByLabelText('Option 1')).toBeInTheDocument();
    expect(screen.getByLabelText('Option 2')).toBeInTheDocument();
  });

  test('does not render options when isExpanded is false', () => {
    render(<FilterSection {...defaultProps} isExpanded={false} />);
    expect(screen.queryByLabelText('Option 1')).not.toBeInTheDocument();
    expect(screen.queryByLabelText('Option 2')).not.toBeInTheDocument();
  });

  test('checkboxes reflect selectedOptions prop', () => {
    const { rerender } = render(<FilterSection {...defaultProps} isExpanded={true} selectedOptions={['testCategory:opt1']} />);
    const checkbox1 = screen.getByLabelText('Option 1');
    const checkbox2 = screen.getByLabelText('Option 2');
    expect(checkbox1).toBeChecked();
    expect(checkbox2).not.toBeChecked();

    // Use rerender to update the props of the already rendered component
    rerender(<FilterSection {...defaultProps} isExpanded={true} selectedOptions={['testCategory:opt2']} />);
    // screen.getByLabelText will now query the updated component
    const checkbox1Updated = screen.getByLabelText('Option 1');
    const checkbox2Updated = screen.getByLabelText('Option 2');
    expect(checkbox1Updated).not.toBeChecked();
    expect(checkbox2Updated).toBeChecked();
  });

  test('calls onFilterChange with correct arguments when a checkbox is clicked', () => {
    render(<FilterSection {...defaultProps} isExpanded={true} />);
    const checkbox1 = screen.getByLabelText('Option 1');
    fireEvent.click(checkbox1);
    expect(mockOnFilterChange).toHaveBeenCalledWith('testCategory', 'opt1');

    const checkbox2 = screen.getByLabelText('Option 2');
    fireEvent.click(checkbox2);
    expect(mockOnFilterChange).toHaveBeenCalledWith('testCategory', 'opt2');
  });
});
