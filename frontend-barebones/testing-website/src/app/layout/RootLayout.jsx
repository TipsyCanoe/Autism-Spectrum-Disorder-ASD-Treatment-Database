import { Outlet } from "react-router-dom";

const RootLayout = () => {
  return (
    <div>
      <h1></h1>
      <Outlet /> {
        
      }
    </div>
  );
};

export default RootLayout;
