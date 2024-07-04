import fabulousLogo from '../public/fabulous.png'
import './App.css'

function App() {

  return (
    <>
      <div>
        <img src={fabulousLogo} className="logo fabulous" alt="Fabulous logo" />
      </div>
      <div>
        <p>Coming to you soon</p>
      </div>
      <p className="read-the-docs">
        Brought to you by a Fabulous team
      </p>
    </>
  )
}

export default App