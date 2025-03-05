import { createBrowserRouter } from "react-router-dom";
import HomePage from "../../routes/HomePage";
import NotFoundPage from "../../routes/NotFoundPage";
import SearchPage from "../../routes/SearchPage";
import SearchResults from "../../routes/SearchResults";
import FAQ from "../../routes/FAQ";
import About from "../../routes/About";
import RootLayout from "../layout/RootLayout";

const router = createBrowserRouter([
  {
    path: "/",
    element: <RootLayout />,
    errorElement: <NotFoundPage />,
    children: [
      { index: true, element: <HomePage /> },
      {
        path: "search",
        element: <SearchPage />,
        children: [{ path: "results", element: <SearchResults /> }],
      },
      { path: "FAQ", element: <FAQ /> },
      { path: "about", element: <About /> },
    ],
  },
]);

export default router;
