import { Grid } from "@mui/joy";

interface OverviewProps {
  identifier: string;
  chainType: string;
  species: string;
  isotypeAndSubclass: string;
}

export default function Overview({
  identifier,
  chainType,
  species,
  isotypeAndSubclass,
}: OverviewProps) {
  return (
    <>
      <Grid
        xs={12}
        sx={{
          border: "1px solid #ccc",
        }}
      >
        Identifier: {identifier}
      </Grid>
      <Grid
        xs={2}
        sx={{
          border: "1px solid #ccc",
        }}
      >
        Chain: {chainType}
      </Grid>
      <Grid
        xs={2}
        sx={{
          border: "1px solid #ccc",
        }}
      >
        {" "}
        Species: {species}
      </Grid>
      <Grid
        xs={2}
        sx={{
          border: "1px solid #ccc",
        }}
      >
        Isotope & Subclass: {isotypeAndSubclass}
      </Grid>
      <Grid xs={6}></Grid>
    </>
  );
}
