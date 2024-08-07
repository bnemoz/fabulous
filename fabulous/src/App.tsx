import { Stack } from "@mui/joy";
import "./App.css";
import Logo from "./components/Logo";
import SearchBox from "./components/Searchbox";

export default function Fabulous() {
  return (
    <>
      <Stack direction="column" alignItems="center">
        <Logo></Logo>
        <SearchBox></SearchBox>
      </Stack>
      <p className="footer">Brought to you by a Fabulous team | 2024</p>
    </>
  );
}
