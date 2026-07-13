mkdir -p dl
cd dl
base="https://suitesparse-collection-website.herokuapp.com/MM/HB"
for m in nos3 nos5 nos6 nos7 gr_30_30 jagmesh1 dwt_918 dwt_992 lshp_577 lshp_778 lshp1009 662_bus 685_bus sherman1; do
  if [ ! -f "$m.mtx" ]; then
    curl -sL --max-time 60 "$base/$m.tar.gz" -o "$m.tar.gz" && tar xzf "$m.tar.gz" 2>/dev/null && find . -name "$m.mtx" -exec mv {} . \; 2>/dev/null
  fi
done
ls -la *.mtx 2>/dev/null
