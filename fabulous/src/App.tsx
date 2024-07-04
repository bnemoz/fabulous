import { useState } from 'react'
import fabulousLogo from '../public/fabulous.png'
import * as React from "react";
import {Image, StyleSheet, View, Text} from "react-native";
import { Padding, Color, FontSize, FontFamily, Border, StyleVariable } from "../GlobalStyles";
import './App.css'

function App() {
  const [count, setCount] = useState(0)

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