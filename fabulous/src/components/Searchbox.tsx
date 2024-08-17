import { Stack } from "@mui/joy";
import Select from "@mui/joy/Select";
import Option from "@mui/joy/Option";
import Button from "@mui/joy/Button";
import Textarea from "@mui/joy/Textarea";
import Link from "@mui/joy/Link";
import axios from "axios";
import { useState } from "react";
import SearchResults from "./searchResults/SearchResults";

export default function SearchBox() {
  const [abSequence, setAbSequence] = useState("");
  const [abDict, setAbDict] = useState(undefined);

  async function search() {
    const response = await axios.get(
      "https://api.fabulous-abs.com/id/?sequence=" + abSequence,
    );
    setAbDict(response.data);
  }

  function populateDummySequence() {
    setAbSequence(
      "GCTGGGTTTTCCTTGTTGCTATTCTCGAGGGTGTCCAGTGTGAGGGCCAGCTGCTGGAGTCTGGAGGAGGACTTGCCCAGCCTGGAGGCTCTCTTAGGCTTAGCTGCACCGCTTCCGGCTTCACCTTCAGCAAGAACGCCATGAACTGGGTGAGGCAGGCCCCTGGAAAAAGGCTGGAGTGGGTGGCCGGCATCATTGGAAACGGCAGCGACACCTACTACGCCGACTCTGTGAAGGGCAGGTTCACCATCAGCAGGGACAACAGCAAGAATACCGTGAGCCTGCAGATGAACTCTCTGAGGGCTGAGGATTCTGCTATTTATTACTGTGCTAAAGACAGACATCCTTGGAGATGGCTGCAGCTTTTTGATAGCTGGGGCCAAGGCACCCTGGTTACAGTTTCTTCTGCTAGCACCAAGGGCCCATCGGTCTTCC",
    );
  }

  return (
    <>
      <Stack direction="row" spacing={2} width={800}>
        <Select defaultValue="human" name="species">
          <Option value="human">Human</Option>
          <Option value="mouse">Mouse</Option>
          <Option value="monkey">Monkey</Option>
          <Option value="humanized">Humanized</Option>
        </Select>

        <Textarea
          name="sequence"
          value={abSequence}
          onChange={(e) => {
            e.stopPropagation();
            setAbSequence(e.target.value);
          }}
          onKeyDown={(e) => {
            if (e.key === "Enter") {
              e.preventDefault();
              search();
            }
          }}
          placeholder="Add one or multiple NT or AA sequence(s) in raw or FASTA format"
          required
          sx={{ width: 600 }}
        />

        <Button color="primary" onClick={search}>
          Search
        </Button>
      </Stack>

      <p className="hint">
        <Link component="a" onClick={populateDummySequence}>
          Don't have any antibody? Click here to load a dummy sequence
        </Link>
      </p>
      <br></br>
      <SearchResults abDict={abDict}></SearchResults>
    </>
  );
}
