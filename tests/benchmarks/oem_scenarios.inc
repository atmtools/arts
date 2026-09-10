// Included after the native OEM regression fixture by run_oem_comparison.py.
CovarianceMatrix scenario_covariance(Index n, std::string_view kind, bool mixed) {
  CovarianceMatrix result;
  const bool sparse = kind == "sparse_inverse" or kind == "diagonal";
  const auto add = [&](Index off, Index count, Index id, bool correlated) {
    Matrix block(count,count,0.);
    for(Index i=0;i<count;++i) for(Index j=0;j<count;++j)
      block[i,j] = (i==j ? 1.+Numeric(i%7)/7. : 0.) +
        (correlated ? 0.15*std::exp(-std::abs(Numeric(i-j))/5.) : 0.);
    if(sparse) {
      Sparse s(count,count);
      for(Index i=0;i<count;++i) for(Index j=0;j<count;++j)
        if(block[i,j]!=0) s.rw(i,j)=block[i,j];
      result.add_correlation({Range(off,count),Range(off,count),{id,id},std::move(s)});
    } else result.add_correlation({Range(off,count),Range(off,count),{id,id},std::move(block)});
  };
  if(not mixed) add(0,n,0,not sparse);
  else {
    const Index a=n/2, b=(n-a)/2, c=n-a-b;
    add(0,a,0,false); add(a,b,1,false); add(a+b,c,2,false);
    Matrix cross(b,c,0.);
    for(Index i=0;i<std::min(b,c);++i) cross[i,i]=0.2;
    if(sparse) {
      Sparse s(b,c);
      for(Index i=0;i<std::min(b,c);++i) s.rw(i,i)=0.2;
      result.add_correlation({Range(a,b),Range(a+b,c),{1,2},std::move(s)});
    } else result.add_correlation({Range(a,b),Range(a+b,c),{1,2},std::move(cross)});
  }
  return result;
}
void scenario_benchmark() {
  using Clock=std::chrono::steady_clock;
  std::cout << "scenario,measurements,states,cache_case,mixed,method,gain,phase,samples,calls_per_sample,median_ms,min_ms,max_ms,status,iterations,cost,sa_inverse_blocks,se_inverse_blocks,state\n";
  Index scenario=0;
  for(const auto [m,n] : {std::pair<Index,Index>{500,100},{50,50},{100,500},{50,10},{5,5},{10,50}})
  for(const auto kind : {"dense_inverse","sparse_inverse","diagonal","cholesky"})
  for(const bool mixed : {false,true}) {
    ++scenario;
    Matrix k(m,n);
    for(Index i=0;i<m;++i) for(Index j=0;j<n;++j)
      k[i,j]=(std::sin((Numeric(i)+1)*(Numeric(j)+1)*0.13)+(i==j ? 2.:0.))/std::sqrt(Numeric(n));
    const auto create=[&] {
      auto owner=std::make_unique<Retrieval>();
      auto& r=*owner;
      r.xa=Vector(n,0.); r.y=Vector(m,0.);
      for(Index i=0;i<m;++i) r.y[i]=std::cos(Numeric(i)*0.17);
      r.sa=scenario_covariance(n,kind,mixed);
      r.se=scenario_covariance(m,kind,mixed);
      r.set_target_size(n);
      r.forward=[&k](const Vector& x, Vector& y, Matrix& jac, bool with_jac) {
        y.resize(k.nrows()); mult(y,k,x);
        if(with_jac) jac=k; else jac.resize(0,0);
      };
      return owner;
    };
    // Retrieval's agenda callback captures this: construct in its final location.
    for(const auto method : {"li","li_m","li_cg","li_cg_m","gn","gn_m","gn_cg","gn_cg_m","lm","lm_cg"})
    for(const bool gain : {false,true})
    for(const bool cold : {true,false}) {
      auto owner=create();
      auto& r=*owner;
      r.clear_matrices=gain ? 0:1;
      const auto sa=r.sa, se=r.se;
      const bool explicit_inverse=std::string_view(kind).ends_with("inverse");
      const auto reset=[&] {
        r.x=Vector{}; r.yf=Vector{}; r.jac=Matrix{}; r.evaluations.clear();
        if(cold) {r.sa=sa;r.se=se;}
      };
      const auto execute=[&] {
        if(explicit_inverse and cold) {r.sa.compute_inverse();r.se.compute_inverse();}
        r.run(method);
      };
      if(explicit_inverse and not cold) {r.sa.compute_inverse();r.se.compute_inverse();}
      reset(); execute();
      const Index count = m<=50 and n<=50 ? 8:1;
      std::vector<double> times;
      for(int sample=0;sample<5;++sample) {
        double ms=0;
        for(Index repeat=0;repeat<count;++repeat) {
          reset(); const auto start=Clock::now();execute();
          ms+=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
          require(r.errors.empty() and r.diagnostics[0]<=1,"OEM scenario failed");
          for(const auto x:r.x) require(std::isfinite(x),"Nonfinite retrieved state");
        }
        times.push_back(ms/Numeric(count));
      }
      std::ranges::sort(times);
      std::cout<<scenario<<','<<m<<','<<n<<','<<kind<<','<<mixed<<','<<method<<','<<gain<<','<<(cold?"cold":"warm")
        <<",5,"<<count<<','<<times[2]<<','<<times.front()<<','<<times.back()<<','<<r.diagnostics[0]<<','<<r.diagnostics[4]<<','<<r.diagnostics[2]
        <<','<<r.sa.ninvdiagblocks()<<','<<r.se.ninvdiagblocks()<<',';
      for(const auto x:r.x) std::cout<<std::setprecision(17)<<x<<';';
      std::cout<<'\n'<<std::flush;
    }
  }
}
