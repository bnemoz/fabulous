import { Stack } from "@mui/joy";
import "./App.css";
import Logo from "./components/Logo";
import SearchBox from "./components/Searchbox";
import About from "./components/About";
import HowItWorks from "./components/HowItWorks";

export default function Fabulous() {
  return (
    <>
      <Stack direction="column" alignItems="center">
        <Logo></Logo>
        <SearchBox></SearchBox>
      </Stack>
      <p className="footer">
        <p><About></About> - <HowItWorks></HowItWorks></p>
        <p>Brought to you by a Fabulous team | 2024</p>
      </p>
    </>
  );
}
