import Grid from "@mui/joy/Grid";
import fabulousLogo from "/fabulous.png";

export default function Logo() {
  return (
    <Grid>
      <a href="https://fabulous.nemoz.me">
        <img src={fabulousLogo} className="logo fabulous" alt="Fabulous logo" />
      </a>
    </Grid>
  );
}
