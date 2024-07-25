import { useState, useEffect } from 'react';
import Box from '@mui/joy/Box';
import Button from '@mui/joy/Button';
import Textarea from '@mui/joy/Textarea';
import fabulousLogo from '/fabulous.png';
import Link from '@mui/joy/Link';
import './App.css';
import { Stack } from '@mui/joy';
import Select from '@mui/joy/Select';
import Option from '@mui/joy/Option';
import axios from 'axios';

function App() {

  const [AbSequence, setAbSequence] = useState('');
  const [searchPerformed, setSearchPerformed] = useState(false);
  const [Species, setSpecies] = useState('human');
  const [antibodydict, setAbDict] = useState([]);

  const handleLinkClick = () => {
    setAbSequence('>Test_Antibody\nEGQLLESGGGLAQPGGSLRLSCTASGFTFSKNAMNWVRQAPGKRLEWVAGIIGNGSDTYYADSVKGRFTISRDNSKNTVSLQMNSLRAEDSAIYYCAKDRHPWRWLQLFDSWGQGTLVTVSS');
  };

  const handleSearch = () => { if (AbSequence !== '') {
    console.log('Submitted sequence and species:', AbSequence, Species);
    setSearchPerformed(true);
  }};

  const fetchAPI = async () => { if (searchPerformed) if (AbSequence) {
    const response = await axios.get("http://127.0.0.1:8282/api/antibody/?sequence="+AbSequence);
    setAbDict(response.data);
  }}

  useEffect(() => {
    fetchAPI();
  })

  if (searchPerformed) {

    return (
      <>
        <div>
          <a href="https://fabulous.nemoz.me">
            <img src={fabulousLogo} className="logo fabulous" alt="Fabulous logo"/>
          </a>
          <p>
            <Button onClick={() => {setSearchPerformed(false), setAbSequence(''), setSpecies('human')}}>Home</Button>
          </p>
        </div>

        <Box sx={{ width: 1200, margin: '20px', padding: '20px', border: '1px solid #ccc' }}>
          {AbSequence}
        </Box>

        <Box sx={{height: 50, width: 120, margin: '20px', padding: '20px', border: '1px solid #ccc' }}>Overview (Antibody ID) 
          Species: {Species}
        </Box>

        <Box sx={{height: 200, width: 120, margin: '20px', padding: '20px', border: '1px solid #ccc' }}>Genes 
          <p>IGH: {}</p>
          <p>IGD: {}</p>
          <p>IGV: {}</p>
          <p>SHM: {}</p>
        </Box>

        <Box sx={{height: 200, width: 120, margin: '20px', padding: '20px', border: '1px solid #ccc' }}>Ab dictionary
          {antibodydict}
          {/* {abdict.map((item, index) => (
            <p key={index}>{item}</p>
          ))} */}
        </Box>
      </>
    );
  }


  return (
    <>
      <div>
        <a href="https://fabulous.nemoz.me">
          <img src={fabulousLogo} className="logo fabulous" alt="Fabulous logo"/>
        </a>
      </div>

      <div>
        <Stack direction={'row'} gap={2} width={800}>

          <Select defaultValue="human" name='species'>
            <Option value="human">Human</Option>
            <Option value="mouse">Mouse</Option>
            <Option value="monkey">Monkey</Option>
            <Option value="humanized">Humanized</Option>
          </Select>

          <Textarea name='sequence'
            value={AbSequence}
            onChange={(e) => setAbSequence(e.target.value)}
            placeholder="Add one or multiple NT or AA sequence(s) in raw or FASTA format"
            required
            sx={{width: 600}}
          />

          <Button color="primary" onClick={handleSearch}>
            Search
          </Button>

        </Stack>

        <p className='hint'>
        <Link component="a" onClick={handleLinkClick}>Don't have any antibody? Click here to load a dummy sequence</Link>
          </p>
        
      </div>

      <Box sx={{height: 50}}></Box>

      <p className="footer">
        Brought to you by a Fabulous team | 2024
      </p>
    </>
  )
}




export default App