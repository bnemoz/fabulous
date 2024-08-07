import { Grid, Stack } from "@mui/joy";
import Overview from "./Overview";
import GeneUsage from "./GeneUsage";
import FullLengthSequences from "./FullLengthSequences";
import RegionsAndSequences from "./RegionsAndSequences";

interface AbDict {
  species: string;
  chain: string;
  isotype: string;
  vdj_aa: string;
  vdj_nt: string;
}

interface SearchResultsProps {
  abSequence: string;
  abDict?: AbDict;
}

export default function SearchResults({
  abSequence,
  abDict,
}: SearchResultsProps) {
  console.debug(abDict);
  if (abDict === null) {
    return;
  }
  return (
    <>
      <Grid container spacing={6}>
        <Overview
          identifier={abSequence}
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
