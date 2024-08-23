import { Grid, Stack } from "@mui/joy";
import Overview from "./Overview";
import GeneUsage from "./GeneUsage";
import FullLengthSequences from "./FullLengthSequences";
import RegionsAndSequences from "./RegionsAndSequences";

interface AbDict {
  oriented_input: string;
  species: string;
  chain: string;
  isotype: string;
  vdj_aa: string;
  vdj_nt: string;
  seq_id: string;
}

export default function SearchResults({ abDict }: { abDict?: AbDict }) {
  if (abDict == undefined) {
    return;
  }
  return (
    <>
      <Grid container spacing={6}>
        <Overview
          identifier={abDict.seq_id}
          chainType={abDict.chain}
          species={abDict.species}
          isotypeAndSubclass={abDict.isotype}
        ></Overview>

        <Grid xs={6}>
          <Stack direction="column" spacing={2}>
            <GeneUsage chain={abDict.chain}></GeneUsage>
            <FullLengthSequences
              aminoAcidVDJ={abDict.vdj_aa}
              nucleotideVDJ={abDict.vdj_nt}
            ></FullLengthSequences>
          </Stack>
        </Grid>

        <Grid xs={6}>
          <Stack direction="column" spacing={2}>
            <RegionsAndSequences></RegionsAndSequences>
          </Stack>
        </Grid>
      </Grid>
    </>
  );
}
