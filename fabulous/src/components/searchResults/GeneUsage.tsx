export default function GeneUsage({ abDict }: { abDict: any }) {
  return (
    <>
      <div style={{width: '100%', height: '100%', background: 'white', borderRadius: 8, border: '1px #DDDDDD solid', flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', display: 'inline-flex'}}>
        <div style={{alignSelf: 'stretch', height: 150, padding: 24, borderRadius: 8, flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', display: 'flex'}}>
            <div style={{alignSelf: 'stretch', justifyContent: 'flex-start', alignItems: 'flex-start', gap: 4, display: 'inline-flex'}}>
                <div style={{flex: '1 1 0', paddingBottom: 24, flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', gap: 8, display: 'inline-flex'}}>
                    <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 24, fontFamily: 'Roboto', fontWeight: '700', lineHeight: 32.02, wordWrap: 'break-word'}}>Gene usage</div>
                </div>
            </div>
            <div style={{alignSelf: 'stretch', paddingLeft: 4, paddingRight: 4, justifyContent: 'flex-start', alignItems: 'flex-start', gap: 32, display: 'inline-flex'}}>
                <div style={{flexDirection: 'column', justifyContent: 'center', alignItems: 'flex-start', gap: 2, display: 'inline-flex'}}>
                    <div style={{color: 'rgba(0, 0, 0, 0.60)', fontSize: 14, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 20.02, letterSpacing: 0.17, wordWrap: 'break-word'}}>V GENE</div>
                    <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 16, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 24, letterSpacing: 0.15, wordWrap: 'break-word'}}>{abDict.v_gene.full}</div>
                </div>
                <div style={{flexDirection: 'column', justifyContent: 'center', alignItems: 'flex-start', gap: 2, display: 'inline-flex'}}>
                    <div style={{color: 'rgba(0, 0, 0, 0.60)', fontSize: 14, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 20.02, letterSpacing: 0.17, wordWrap: 'break-word'}}>D GENE</div>
                    <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 16, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 24, letterSpacing: 0.15, wordWrap: 'break-word'}}>{abDict.d_gene.full}</div>
                </div>
                <div style={{flexDirection: 'column', justifyContent: 'center', alignItems: 'flex-start', gap: 2, display: 'inline-flex'}}>
                    <div style={{color: 'rgba(0, 0, 0, 0.60)', fontSize: 14, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 20.02, letterSpacing: 0.17, wordWrap: 'break-word'}}>J GENE</div>
                    <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 16, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 24, letterSpacing: 0.15, wordWrap: 'break-word'}}>{abDict.j_gene.full}</div>
                </div>
            </div>
        </div>
        <div style={{justifyContent: 'flex-end', alignItems: 'flex-start', display: 'inline-flex'}}>
            <div style={{height: 172, boxShadow: '0px 2px 4px -1px rgba(0, 0, 0, 0.20)', justifyContent: 'flex-end', alignItems: 'center', display: 'flex'}}>
                <div style={{height: 172, opacity: 0, justifyContent: 'flex-end', alignItems: 'center', display: 'flex'}}>
                    <div style={{flex: '1 1 0', paddingTop: 16, paddingBottom: 16, paddingLeft: 20, paddingRight: 32, background: '#616161', borderRadius: 4, overflow: 'hidden', flexDirection: 'column', justifyContent: 'center', alignItems: 'flex-start', gap: 12, display: 'inline-flex'}}>
                        <div style={{width: 16, height: 16, padding: 3.33, justifyContent: 'center', alignItems: 'center', display: 'inline-flex'}}>
                            <div style={{width: 9.33, height: 9.33, background: 'white'}}></div>
                        </div>
                        <div style={{alignSelf: 'stretch', color: 'white', fontSize: 10, fontFamily: 'Roboto', fontWeight: '500', lineHeight: 14, wordWrap: 'break-word'}}>Somatic hypermutation rates are calculated as the distance to V-gene or V- and J-genes germlines, with output in nucleotide or amino-acids. These calculations are based on the initial assignment by AbStar. </div>
                    </div>
                    <div style={{width: 12, height: 6, transform: 'rotate(90deg)', transformOrigin: '0 0', background: 'rgba(97, 97, 97, 0.90)'}}></div>
                </div>
            </div>
            <div style={{padding: 8, borderRadius: 100, overflow: 'hidden', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', display: 'inline-flex'}}>
                <div style={{justifyContent: 'flex-start', alignItems: 'flex-start', display: 'inline-flex'}}>
                    <div style={{width: 24, height: 24, position: 'relative'}}>
                        <div style={{width: 19.50, height: 19.50, left: 2.25, top: 2.25, position: 'absolute', background: 'rgba(0, 0, 0, 0.56)'}}></div>
                    </div>
                </div>
            </div>
        </div>
        <div style={{alignSelf: 'stretch', height: 245, padding: 24, flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', display: 'flex'}}>
            <div style={{alignSelf: 'stretch', justifyContent: 'flex-start', alignItems: 'flex-start', gap: 4, display: 'inline-flex'}}>
                <div style={{flex: '1 1 0', paddingBottom: 24, flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', gap: 8, display: 'inline-flex'}}>
                    <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 24, fontFamily: 'Roboto', fontWeight: '700', lineHeight: 32.02, wordWrap: 'break-word'}}>Somatic Hypermutation</div>
                    <div style={{alignSelf: 'stretch', color: 'rgba(0, 0, 0, 0.60)', fontSize: 16, fontFamily: 'Inter', fontWeight: '400', lineHeight: 19.20, wordWrap: 'break-word'}}>Distance to germline</div>
                </div>
            </div>
            <div style={{alignSelf: 'stretch', height: 114, flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', gap: 2, display: 'flex'}}>
                <div style={{alignSelf: 'stretch', height: 42, flexDirection: 'column', justifyContent: 'flex-start', alignItems: 'flex-start', display: 'flex'}}>
                    <div style={{alignSelf: 'stretch', justifyContent: 'flex-start', alignItems: 'flex-start', display: 'inline-flex'}}>
                        <div style={{flex: '1 1 0', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', display: 'inline-flex'}}>
                            <div style={{alignSelf: 'stretch', paddingTop: 9, paddingBottom: 9, paddingLeft: 25, paddingRight: 16, justifyContent: 'center', alignItems: 'center', gap: 8, display: 'inline-flex'}}>
                                <div style={{color: '#9C27B0', fontSize: 14, fontFamily: 'Roboto', fontWeight: '500', textTransform: 'uppercase', lineHeight: 24, letterSpacing: 0.40, wordWrap: 'break-word'}}>V gene  </div>
                            </div>
                            <div style={{width: 286, height: 0, background: '#9C27B0', border: '2px #9C27B0 solid'}}></div>
                        </div>
                        <div style={{flex: '1 1 0', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', display: 'inline-flex'}}>
                            <div style={{paddingLeft: 16, paddingRight: 16, paddingTop: 9, paddingBottom: 9, justifyContent: 'center', alignItems: 'center', gap: 8, display: 'inline-flex'}}>
                                <div style={{color: 'rgba(0, 0, 0, 0.60)', fontSize: 14, fontFamily: 'Roboto', fontWeight: '500', textTransform: 'uppercase', lineHeight: 24, letterSpacing: 0.40, wordWrap: 'break-word'}}>v and j gene</div>
                            </div>
                        </div>
                    </div>
                </div>
                <div style={{alignSelf: 'stretch', justifyContent: 'space-between', alignItems: 'flex-start', display: 'inline-flex'}}>
                    <div style={{flex: '1 1 0', padding: 12, borderRadius: 4, flexDirection: 'column', justifyContent: 'center', alignItems: 'center', gap: 4, display: 'inline-flex'}}>
                        <div style={{border: '1px solid', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', gap: 2, display: 'flex'}}>
                            <div style={{color: 'rgba(0, 0, 0, 0.60)', fontSize: 14, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 20.02, letterSpacing: 0.17, wordWrap: 'break-word'}}>nucleotides</div>
                            <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 16, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 24, letterSpacing: 0.15, wordWrap: 'break-word'}}>{abDict.nt_identity.V}%</div>
                        </div>
                    </div>
                    <div style={{width: 65, height: 0, transform: 'rotate(90deg)', transformOrigin: '0 0', border: '1px #DDDDDD solid'}}></div>
                    <div style={{flex: '1 1 0', padding: 12, borderRadius: 4, flexDirection: 'column', justifyContent: 'center', alignItems: 'center', gap: 4, display: 'inline-flex'}}>
                        <div style={{border: '1px solid', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', gap: 2, display: 'flex'}}>
                            <div style={{color: 'rgba(0, 0, 0, 0.60)', fontSize: 14, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 20.02, letterSpacing: 0.17, wordWrap: 'break-word'}}>amino-acids</div>
                            <div style={{color: 'rgba(0, 0, 0, 0.87)', fontSize: 16, fontFamily: 'Roboto', fontWeight: '400', lineHeight: 24, letterSpacing: 0.15, wordWrap: 'break-word'}}>{abDict.aa_identity.V}%</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        <div style={{justifyContent: 'flex-end', alignItems: 'flex-start', display: 'inline-flex'}}>
            <div style={{height: 172, boxShadow: '0px 2px 4px -1px rgba(0, 0, 0, 0.20)', justifyContent: 'flex-end', alignItems: 'center', display: 'flex'}}>
                <div style={{height: 172, opacity: 0, justifyContent: 'flex-end', alignItems: 'center', display: 'flex'}}>
                    <div style={{flex: '1 1 0', paddingTop: 16, paddingBottom: 16, paddingLeft: 20, paddingRight: 32, background: '#616161', borderRadius: 4, overflow: 'hidden', flexDirection: 'column', justifyContent: 'center', alignItems: 'flex-start', gap: 12, display: 'inline-flex'}}>
                        <div style={{width: 16, height: 16, padding: 3.33, justifyContent: 'center', alignItems: 'center', display: 'inline-flex'}}>
                            <div style={{width: 9.33, height: 9.33, background: 'white'}}></div>
                        </div>
                        <div style={{alignSelf: 'stretch', color: 'white', fontSize: 10, fontFamily: 'Roboto', fontWeight: '500', lineHeight: 14, wordWrap: 'break-word'}}>Somatic hypermutation rates are calculated as the distance to V-gene or V- and J-genes germlines, with output in nucleotide or amino-acids. These calculations are based on the initial assignment by AbStar. </div>
                    </div>
                    <div style={{width: 12, height: 6, transform: 'rotate(90deg)', transformOrigin: '0 0', background: 'rgba(97, 97, 97, 0.90)'}}></div>
                </div>
            </div>
            <div style={{padding: 8, borderRadius: 100, overflow: 'hidden', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', display: 'inline-flex'}}>
                <div style={{justifyContent: 'flex-start', alignItems: 'flex-start', display: 'inline-flex'}}>
                    <div style={{width: 24, height: 24, position: 'relative'}}>
                        <div style={{width: 19.50, height: 19.50, left: 2.25, top: 2.25, position: 'absolute', background: 'rgba(0, 0, 0, 0.56)'}}></div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    </>
  );
}
