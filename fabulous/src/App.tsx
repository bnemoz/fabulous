import { useState } from 'react';
import Box from '@mui/joy/Box';
import Button from '@mui/joy/Button';
import Textarea from '@mui/joy/Textarea';
import fabulousLogo from '/fabulous.png';
import Link from '@mui/joy/Link';
import './App.css';
import { Stack } from '@mui/joy';
import Select from '@mui/joy/Select';
import Option from '@mui/joy/Option';

function App() {

  const [AbSequence, setAbSequence] = useState('');

  const handleLinkClick = () => {
    setAbSequence('>Test_Antibody\nEGQLLESGGGLAQPGGSLRLSCTASGFTFSKNAMNWVRQAPGKRLEWVAGIIGNGSDTYYADSVKGRFTISRDNSKNTVSLQMNSLRAEDSAIYYCAKDRHPWRWLQLFDSWGQGTLVTVSS');
  };

  return (
    <>
      <div>
        <img src={fabulousLogo} className="logo fabulous" alt="Fabulous logo"/>
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

          <Button color="primary">
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