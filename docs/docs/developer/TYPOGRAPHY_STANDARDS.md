# Typography Standards

This document outlines the standardized font sizes and typography utilities for consistent styling across the frontend application.

## Custom CSS Classes (Recommended)

These classes are defined in `frontend/testing-website/src/index.css` and provide consistent styling across all pages:

### Headings
- **`.page-heading`** - Main page titles (e.g., "About the Autism Resources Database")
  - Size: `text-3xl` (30px)
  - Weight: Bold
  - Usage: One per page, centered

- **`.section-heading`** - Section headings (e.g., "Our Mission", "What We Offer")
  - Size: `text-2xl` (24px)
  - Weight: Semibold
  - Usage: Multiple per page

- **`.card-heading`** - Card/component titles
  - Size: `text-xl` (20px)
  - Weight: Semibold
  - Usage: Feature cards, resource cards, etc.

### Body Text
- **`.body-text`** - Standard body text for main content
  - Size: `text-lg` (18px)
  - Usage: Paragraphs, descriptions

- **`.body-text-sm`** - Smaller body text for secondary content
  - Size: `text-base` (16px)
  - Usage: Captions, metadata, helper text

### Buttons
- **`.btn-standard`** - Base button styling
  - Padding: `px-6 py-2.5`
  - Size: `text-base` (16px)
  - Border radius: rounded-lg

- **`.btn-primary`** - Primary action buttons
  - Extends `.btn-standard`
  - Background: navbar-blue
  - Hover: link-hover-blue

## Direct Tailwind Classes (Alternative)

If you prefer not to use the custom classes, here's the standard sizing:

### Headings
```css
h1: text-3xl font-bold
h2: text-2xl font-semibold
h3: text-xl font-semibold
```

### Body Text
```
Paragraphs: text-lg
Small text: text-base
Captions: text-sm
```

### Buttons
```
Standard: px-6 py-2.5 text-base
Large: px-8 py-3 text-lg
```

## Current Implementation

### Pages Using Standards:
- ✅ **About.jsx** - Uses custom CSS classes
- ✅ **FAQ.jsx** - Uses custom CSS classes
- ✅ **FAQItem.jsx** - Uses custom CSS classes
- ✅ **HomePage.jsx** - Uses custom CSS classes

### Usage Example:

```jsx
// Instead of:
<h1 className="text-3xl font-bold text-center mb-8">
  About the Autism Resources Database
</h1>

// Use:
<h1 className="page-heading">
  About the Autism Resources Database
</h1>
```

## Benefits of Standardization

1. **Consistency** - All pages look cohesive
2. **Easy Updates** - Change font sizes in one place (index.css)
3. **Responsive Design** - Add responsive variants to custom classes
4. **Maintainability** - Developers know which class to use
5. **Reduced Code** - Less repetitive Tailwind classes

## File Locations

- **CSS Classes**: `frontend/testing-website/src/index.css`
- **Implementation**: All page components in `frontend/testing-website/src/routes/`

## Future Improvements

Consider adding responsive variants for better mobile experience:

```css
.page-heading {
  @apply text-2xl md:text-3xl lg:text-4xl font-bold text-center mb-8;
}
```

This would make headings progressively larger on tablets and desktops automatically.

## Testing

All typography changes are covered by frontend tests:
- `HomePage.test.js` - Validates homepage styling
- `App.test.js` - Validates navigation styling
- Component tests validate proper class usage

Run tests with:
```bash
cd frontend/testing-website
npm test
```
