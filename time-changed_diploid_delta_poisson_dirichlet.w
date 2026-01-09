%$ echo "const unsigned int SAMPLE_SIZE = 50;" > headerfile.hpp
%$ echo "const double RHO = 100.0;" >> headerfile.hpp
%$ echo "const double ALPHA = 0.01;" >> headerfile.hpp
%$ echo "const double KAPPA = 2;" >> headerfile.hpp
%$ echo "const double CEPS = 100;" >> headerfile.hpp
%$ echo "const double MINF = (2. + (1 + pow(2., 1-KAPPA)/(KAPPA-1)))/2. ;"  >>headerfile.hpp
%$ echo "const double BKAPPA = pow(2,KAPPA)*(KAPPA - 2)*(KAPPA - 1);" >> headerfile.hpp
%$ echo "const double AKAPPA = ((KAPPA+2) + pow(KAPPA,2))/2.;" >> headerfile.hpp
%$ echo "const double CKAPPA = ((KAPPA > 2 ? AKAPPA / BKAPPA : 1.)*2.) / pow(MINF,2.) ;" >> headerfile.hpp
%$ echo "const int EXPERIMENTS = 1e5;" >> headerfile.hpp
%$ NAFN=time-changed_diploid_delta_poisson_dirichlet
%$ ctangle $NAFN
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 484843-30383820202 -n1) > logitresout
%$ seq 49 | awk '{n=50;print log($1/n) - log(1 - ($1/n))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ NAFN=time-changed_diploid_delta_poisson_dirichlet
%$ cweave $NAFN
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ emacs --version | head -n1 > innleggemacs
%$ g++ --version | head -n1 > innleggcpp
%$ lualatex --version | head -n2 > innlegglualatex
%$ cweave  --version | head -n1 > innleggcweave
%$ ctangle  --version | head -n1  > innleggctangle
%$ uname  --kernel-release -o  > innleggop
%$ bash --version | head -n1 > innleggbash
%$ sed -i 's/x86/x86\\/g' innleggbash
%$ gsl-config --version > innlegggsl
%$ spix --version > innleggspix
%$ lualatex $NAFN
%$ lualatex $NAFN
\documentclass[a4paper,10pt,final]{cweb}
\usepackage{fancyvrb,  graphicx, url}
\usepackage[svgnames]{xcolor}
\usepackage[inputlexerlinenos=true]{minted}
\setminted[cpp]{breaklines}
\usepackage{piton}
\PitonOptions{width=min,line-numbers,break-lines,indent-broken-lines,background-color=gray!5}
%% a4 paper size 210 x 297 millimeters
%% for xelatex
%%\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
%%\usepackage[all]{xy}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage{natbib}
\usepackage{abstract}
\renewcommand{\abstractname}{}
%\usepackage[backend=biber,style=authoryear-icomp,
%    sortlocale=de_DE,
%    natbib=true,
%    url=false, 
%    doi=true,
%    eprint=false]{biblatex}
%    \addbibresource{refs.bib}
    \usepackage[]{hyperref}
 %%   \hypersetup{colorlinks=true}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{bbm}
\usepackage{orcidlink}
\usepackage[en-GB,showdow]{datetime2}
\DTMlangsetup[en-GB]{ord=raise,monthyearsep={,\space}}
\usepackage[rightcaption]{sidecap}
\sidecaptionvpos{figure}{m}
\usepackage{siunitx}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=newest}
\usepgfplotslibrary{units}
\sisetup{
  round-mode          = places,
  round-precision     = 3,
}
\setstretch{1}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\IN}{\ensuremath{\mathds{N}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\textbf{\small #1}}}
\newcommand{\s}{{\textcolor{blue}{\ensuremath\bullet}}}
\newcommand{\sa}[1]{\mathrel{\reflectbox{\rotatebox[origin=c]{#1}{\scalebox{0.5}{$\uparrow$}}}}}
\newcommand{\x}{{\circ}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  {\textsf{\large \textbf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\large \@@author}
\end{center}
\egroup
}
\makeatother
\title{Gene genealogies in diploid populations evolving according to
sweepstakes reproduction \\ --- approximating $\EE{R_i(n)}$ for the
 time-changed $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent}
 \author{Bjarki Eldon
 \footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{ compiled @@
{\DTMcurrenttime} on  {\today}  \\
\input{innleggctangle} \\
 \input{innleggcpp} \\ kernel  \input{innleggop} \\  \input{innleggbash} \\
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  {\LaTeX}
\input{innlegglualatex}  \\  written using  \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391}}


\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

\begin{abstract}
Let $\# A$ denote the number of elements in a finite set $A$.  For a
given coalescent $\set{\xi^{n}}$ write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :
\#\xi = i }dt$ and $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $
where $\tau(n) \equiv \inf \set{t \ge 0 : \#\xi^{n}(t) = 1}$. Then
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Write
$R_{i}(n) \equiv L_{i}(n)/L(n)$ for $i = 1,2,\ldots, n-1$.  With this
C++ simulation code we approximate $\EE{R_{i}(n)}$ when
$\set{\xi^{n}(G)} \equiv \set{\xi^{n}(G(t)) : t \ge 0 } $ is the
time-changed continuous-time  $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$
coalescent with $0 < \alpha < 1$ and time-change function $G(t)$.
\end{abstract}

\tableofcontents


@* {\bf Copyright}.
\label{sec:copy}


Copyright {\textcopyright} {\the\year} Bjarki Eldon \newline 

This document and any source code it contains is distributed under the
terms of the GNU General Public Licence (version $\ge 3$).  You should
have received a copy of the licence along with this file (see file
COPYING).


The source codes described in this document are free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This document and the code it contains is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU General Public License for more details.  You should have
received a copy of the GNU General Public License along with this file
(see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf compilation and output}. 
\label{sec:comp}

This CWEB \citep{knuth1994cweb} document (the {\tt .w} file) can be
compiled with {\tt cweave} to generate a {\tt .tex} file, and with
{\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c} C++ code
file. 

Use the shell tool {\tt spix} on the  script appearing before  the
preamble (the lines starting with  \%\$); simply 

{\tt spix /path/to/the/sourcefile}

where {\tt sourcefile} is the {\tt .w} file

One   may also  copy the script into a file and run  {\tt parallel}
\cite{tange11:_gnu_paral} :

{\tt parallel --gnu -j1 :::: /path/to/scriptfile}


@* {\bf intro}.
\label{sec:intro}


Here we consider the  continuous-time $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent
corresponding to exponential population growth where the time-change
function $G(t) = \int_{0}^{t}e^{\rho s}ds $ for some fixed $\rho \ge
0$.


One samples  the
$\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent by
sampling group sizes $2 \le k_{1},\ldots,k_{r} \le n$ with
$\sum_{i}k_{i}\le n$, $s = n - \sum_{i}k_{i}$,    according to  
\begin{equation}
\label{eq:1}
\begin{split}
\lambda_{n;k_{1},\ldots,k_{r};s} & = \one{r=1,k_{1}=2} \frac{C_{\kappa}}{C_{\kappa} + c(1-\alpha)} +  \frac{cp_{n;k_{1},\ldots,k_{r};s}}{C_{\kappa} + c(1-\alpha)}\\
C_{\kappa} &  =  \frac{2}{m_{\infty}^{2}}\svigi{\one{\kappa = 2} +  \one{\kappa > 2} \frac{c_{\kappa}}{2^{\kappa}(\kappa - 2)(\kappa - 1) } }  \\
p_{n;k_{1},\ldots,k_{r};s} & =  \frac{\alpha^{r+s-1}\svigi{r+s-1}!}{(n-1)! }\prod_{i=1}^{r}\svigi{k_{i}-1-\alpha}_{k_{i}-1}
\end{split}
\end{equation}
where $2 + \kappa < c_{\kappa} < \kappa^{2}$ when $\kappa > 2$.  
The blocks in each group are then  split among four subgroups
independently and uniformly at random, and  the blocks in the same
subgroup are merged.



Let $\# A$ denote the number of elements in a finite set $A$.  For a
given coalescent $\set{\xi^{n}}$ write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :
\#\xi = i }dt$ and $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $
where $\tau(n) \equiv \inf \set{t \ge 0 : \#\xi^{n}(t) = 1}$. Then
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Write
$R_{i}(n) \equiv L_{i}(n)/L(n)$ for $i = 1,2,\ldots, n-1$.  With this
C++ code we use 
simulations to approximate $\EE{R_{i}(n)}$ when the coalescent is the
time-changed continuous-time
$\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent with
time-change function $G(t) = \int_{0}^{t}e^{\rho s}ds$ for a given
fixed $\rho \ge 0$. 


The code follows in \S~\ref{sec:includes}--\S~\ref{sec:main}, we
conclude in \S~\ref{sec:bib}. Comments within the code in \aths{this
font and color}




@* {\bf code}.
\label{sec:code}



@*1 {\bf includes}.
\label{sec:includes}


the included libraries and header file with the global constants


@<includes@>=@#
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <random>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <assert.h>
#include <float.h>
#include <fenv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <boost/math/special_functions/factorials.hpp>
#include "headerfile.hpp"



@*1 {\bf random number generators}.
\label{sec:rngs}



@<rngs@>=@#
 /* \newline \aths{define a random seed object} */
  std::random_device randomseed;
  /* \newline \aths{ Standard Mersenne twister  random number engine} */
  std::mt19937_64 rng(randomseed());



gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf $e^{x}$ with checks}.
\label{sec:veldi}


compute $e^{x}$ checking for over- and underflow


@<$e^{x}$@>=@#
static double veldi( const double x, const double y )
  {
    feclearexcept(FE_ALL_EXCEPT);
    const double d = pow(x,y);

    return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;   
  }



@*1 {\bf $\lambda_{n;k_{1},\ldots,k_{r};s}$}.
\label{sec:lambdanks}


compute  $\lambda_{n;k_{1},\ldots,k_{r};s}$ \eqref{eq:1}

@<$\lambda_{n;k_{1},\ldots,k_{r};s}$ \eqref{eq:1}@>=@#
static double lambdanks (const double n, const std::vector<unsigned int>& v_k)
{

@#
  assert( v_k.size() > 0);
  @#
  assert(std::all_of( v_k.cbegin(), v_k.cend(), [](const auto k){return k > 1;}));
@#
  
  double d {};
  double k {} ;
  double f {1} ;
  const double r = static_cast<double>(v_k.size());
  @#
  std::unordered_map<unsigned int, unsigned int> counts {} ;
   /* \newline \aths{$(x)_{m} \equiv x(x-1)\cdots (x-m+1)$} */
   auto ff = [](const double x, const unsigned int m){
   return static_cast<double>(boost::math::falling_factorial (x, m));};
  
  for( std::size_t i = 0; i < v_k.size(); ++i){
    f *= ff( static_cast<double>(v_k[i]) - 1. - ALPHA, v_k[i] -1);
    /* \newline \aths{ count occurrence of each  merger size} */
    ++counts[v_k[i]];
    k += static_cast<double>(v_k[i]) ;
    d += lgamma(static_cast<double>(v_k[i] + 1)) ; }

  assert( k < n + 1 ); 
  const double s = n - k;
  
  const double p = static_cast<double>(std::accumulate( counts.cbegin(), counts.cend(), 0, [](double a, const auto &x){ return a +  lgamma( (double)x.second + 1);}));
  
  const double l = ((v_k.size() < 2 ? (v_k[0] < 3 ? 1. : 0) : 0)*CKAPPA) + (CEPS * veldi(ALPHA, r+s-1) * tgamma(r+s)   *f/tgamma(n) ) ;
      /* \newline \aths{\S~\ref{sec:veldi}} */  
  return ( veldi(exp(1), (lgamma( n + 1.) - d)  - lgamma( n - k + 1 ) - p) * l / (CKAPPA + (CEPS*(1-ALPHA))) ) ;
}



@*1 {\bf generate group sizes}.
\label{sec:GenPartitions}


generate group sizes $k_{1},\ldots,k_{r}$ summing to  |myInt| 


@<sizes of groups@>=@#
static double  GenPartitions(const unsigned int m,  const unsigned int myInt,
                   const unsigned int PartitionSize,
                   unsigned int MinVal,
			  unsigned int MaxVal,
			  std::vector< std::pair< double, std::vector<unsigned int> > > & v_l_k,
			  std::vector<double>& lrates_sorting)
{


    /* \newline \aths{ |m| is the given number of blocks;
     the partitions sum to  |myInt|} */ 

  double lrate {} ;
  double sumrates {} ;
  
    std::vector<unsigned int> partition(PartitionSize);
    unsigned int idx_Last = PartitionSize - 1;
    unsigned int idx_Dec = idx_Last;   
    unsigned int idx_Spill = 0;        
    unsigned int idx_SpillPrev;        

    unsigned int LeftRemain = myInt - MaxVal - (idx_Dec - 1)*MinVal;  
    partition[idx_Dec] = MaxVal + 1;    

    

    do {
        unsigned int val_Dec = partition[idx_Dec] - 1;      
        partition[idx_Dec] = val_Dec;                       

        idx_SpillPrev = idx_Spill;          
        idx_Spill = idx_Dec - 1;            

        while (LeftRemain > val_Dec)        
        {
            partition[idx_Spill--] = val_Dec;
            LeftRemain -= val_Dec - MinVal; 
            
        }   

	
        partition[idx_Spill] = LeftRemain;  
        

        const char a = (idx_Spill) ? ~((-3 >> (LeftRemain - MinVal)) << 2) : 11;  
        const char b = (-3 >> (val_Dec - LeftRemain));

        switch (a & b)  
        {
            case 1:
            case 2:
            case 3: idx_Dec = idx_Spill;
                    LeftRemain = 1 + (idx_Spill - idx_Dec + 1)*MinVal; 
                    break;

            case 5: for (++idx_Dec, LeftRemain = (idx_Dec - idx_Spill)*val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= MinVal); idx_Dec++) 
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;

            case 6:
            case 7:
            case 11:idx_Dec = idx_Spill + 1;
                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;


            case 9: for (++idx_Dec, LeftRemain = idx_Dec * val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec + 1)); idx_Dec++)  
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 - (idx_Dec - 1)*MinVal;
                    break;

            case 10:for (LeftRemain += idx_Spill * MinVal + (idx_Dec - idx_Spill)*val_Dec + 1, ++idx_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec - 1)); idx_Dec++)    
                        LeftRemain += partition[idx_Dec];

                    LeftRemain -= (idx_Dec - 1)*MinVal;
                    break;
        }

	
        while (idx_Spill > idx_SpillPrev)   
            partition[--idx_Spill] = MinVal;
            @#
        assert( static_cast<unsigned int>(std::accumulate( partition.begin(), partition.end(), 0)) == myInt);

@#

        lrate =  lambdanks( static_cast<double>(m),  partition) ;

@#
	assert( lrate >= 0) ;
	
	v_l_k.push_back( std::make_pair( lrate, partition) ) ;

	lrates_sorting.push_back( lrate );

	sumrates += lrate ; 
@#	

    } while (idx_Dec <= idx_Last);
@#
    assert( sumrates >= 0) ;
    return sumrates ;
}




@*1 {\bf group sizes up to a given number}.
\label{sec:allmergers_sum_m}

get all partitions summing to a given number

@<sizes up to number@>=@#
static double allmergers_sum_m( const unsigned int n, const unsigned int m, std::vector< std::pair< double, std::vector<unsigned int> > >& v__l_k,  std::vector<double>&  v_lrates_sort )
{
  /* \newline \aths{ |n| is number of blocks ;
      all partitions summing to $|m| \le  |n|$} */
  const std::vector<unsigned int> v__m {m};
  /* \newline \aths{\S~\ref{sec:lambdanks}} */
  double sumr = lambdanks( static_cast<double>(n), v__m ) ;

  v__l_k.push_back( std::make_pair( sumr, v__m ) ) ;
@#
  v_lrates_sort.push_back( sumr );

  if( m > 3){
    for( unsigned int s = 2; s <= m/2; ++s){
      
      assert(m > 2*(s-1) );
      /* \newline \aths{\S~\ref{sec:GenPartitions}} */
      sumr += GenPartitions(n,  m, s, 2, m - (2*(s-1)), v__l_k, v_lrates_sort ); }
  }
  @#
  assert( sumr >= 0) ;

@#
  return sumr ;
}


@*1 {\bf write partitions to files}. 
\label{sec:ratesmergersfile}


write partitions to files {\tt gg\_{n}\_.txt}


@<partitions to files@>=@#
static void ratesmergersfile( const unsigned int n, const std::vector<unsigned int>& v__indx, const std::vector< std::pair< double, std::vector<unsigned int> > > & vlk, const double s,   std::vector< std::vector< double > > & a__cmf )
{
@#
  assert( s > 0);
  double cmf {} ;
  std::ofstream f ;
  @#
  f.open("gg_" + std::to_string(n) + "_.txt", std::ios::app);

@#
  a__cmf[n].clear() ;

@#

  for( const auto &i : v__indx)
    {
      cmf += (vlk[ i ].first) / s ;
      assert( cmf >= 0) ;
     @#
      a__cmf[n].push_back(cmf) ;
      assert((vlk[ i ].second).size() > 0 );
        @#
      for( const auto &x : vlk[ i ].second){
	f << x << ' ' ;}
      f << '\n' ; }
@#

f.close();
@#

  assert(abs(cmf - 1.) < 0.999999);
}


@*1 {\bf partitions when given number of blocks}.
\label{sec:allmergers_when_n_blocks}


@<all partitions given number of blocks@>=@#
static void allmergers_when_n_blocks( const unsigned int n, std::vector<double> & v__lambdan, std::vector< std::vector< double > > & a__cmf )
{
@#
  std::vector< std::pair< double, std::vector<unsigned int> > > vlk {} ;
  std::vector< double > ratetosort {} ;
  ratetosort.clear() ;

  double lambdan {} ;
  vlk.clear() ;
  assert( n > 1);
  for( unsigned int k = 2 ; k  <= n ; ++k){
     /* \newline \aths{ the partition sums to |k|; the number of
     blocks is |n|; \S~\ref{sec:allmergers_sum_m}} */
   lambdan += allmergers_sum_m(n, k, vlk, ratetosort); }
  /* \newline \aths{ record the total rate when |n| blocks;  use for sampling time} */
  assert( lambdan > 0);
  v__lambdan[n] = lambdan ;

@#
  
  std::vector<unsigned int> indx (ratetosort.size());

@#

  std::iota( indx.begin(), indx.end(), 0);

@#

  std::stable_sort( indx.begin(), indx.end(), [&ratetosort](const unsigned int x, const unsigned int y){return ratetosort[x] > ratetosort[y];});

   /* \newline \aths{\S~\ref{sec:ratesmergersfile}} */
  ratesmergersfile(n, indx, vlk,  v__lambdan[n], a__cmf);
}



@*1 {\bf all partitions}.
\label{sec:all_partitions_rates}


generate all partitions and rates

@<partitions generate all@>=@#
static void all_partitions_rates( std::vector<double>& vlmn, std::vector< std::vector<double> > & acmf  )
{
   
  for( unsigned int tmpn = 2; tmpn <= SAMPLE_SIZE; ++tmpn )
    {
      /* \newline \aths{\S~\ref{sec:allmergers_when_n_blocks}} */
      allmergers_when_n_blocks( tmpn, vlmn, acmf );
    }
}




@*1 {\bf read in partitions}.
\label{sec:readmergersizes}


read in partitions  $k_{1},\ldots,k_{r}$ from files {\tt g\_<m>\_.txt}


@<partitions read in@>=@#
static void readmergersizes(const unsigned int n,  const unsigned int j, std::vector<unsigned int> & v__mergers )
{
@#
  std::ifstream f("gg_" + std::to_string(n) + "_.txt") ;

@#

  std::string line {} ;

@#

  v__mergers.clear(); 

@#
  for( unsigned int i = 0; std::getline(f, line) && i < j; ++i ){
    if(i >= j-1 ){
  @#    
      std::stringstream ss(line) ;
      @#
      v__mergers = std::vector<unsigned int>( std::istream_iterator<unsigned int>(ss), {} );
    } }
@#
  assert(v__mergers.size() > 0) ;

@#

  assert( std::all_of(v__mergers.cbegin(), v__mergers.cend(), [](const auto &x){ return x > 1; }) );

@#

  f.close();
}


@*1 {\bf split partitions}.
\label{sec:split_blocks}


\begin{Piton}
void split_blocks(const unsigned int k, std::vector<unsigned int>& split_partition )
{
  auto sample_box = [](void)
  { const double u = gsl_rng_uniform(rngtype);
    return static_cast<unsigned int>(u < 0.25 ? 0 : ( u < .5 ? 1 : ( u < .75 ? 2 : 3))) ;
    } ;

  std::vector<unsigned int> boxes (4);
  for( unsigned int j = 0; j < k; ++j)
    {
      ++boxes[ sample_box()  ];
    }
  assert( static_cast<unsigned int>( std::accumulate(boxes.cbegin(), boxes.cend(), 0))  == k);

  const auto t = static_cast<std::size_t>( std::count_if( boxes.cbegin(), boxes.cend(), [](const auto b){ return b > 1;}));

  if(t > 0){
  std::copy_if(boxes.cbegin(), boxes.cend(), std::back_inserter(split_partition), [](const auto b){return b > 1;});
  }
}
\end{Piton}


@<split partition into boxes@>=@#
void split_blocks(const unsigned int k, std::vector<unsigned int>& split_partition )
{
@#
  auto sample_box = [](void)
  { const double u = gsl_rng_uniform(rngtype);
    return static_cast<unsigned int>(u < 0.25 ? 0 : ( u < .5 ? 1 : ( u < .75 ? 2 : 3))) ;
    } ;

@#
  std::vector<unsigned int> boxes (4);
  @#
  for( unsigned int j = 0; j < k; ++j)
    {
      ++boxes[ sample_box()  ];
    }
@#
  assert( static_cast<unsigned int>( std::accumulate(boxes.cbegin(), boxes.cend(), 0))  == k);

  /* \newline \aths{ count how many boxes have at least 2 blocks} */
  const auto t = static_cast<std::size_t>( std::count_if( boxes.cbegin(), boxes.cend(), [](const auto b){ return b > 1;}));
@#
  if(t > 0){
  std::copy_if(boxes.cbegin(), boxes.cend(), std::back_inserter(split_partition), [](const auto b){return b > 1;});
  }
}


@*1 {\bf split groups in partition}.
\label{sec:split_groups}

 given a partition, return with each  group split into boxes 

\begin{Piton}
void split_groups( const std::vector< unsigned int>& partition, std::vector<unsigned int>& split_partition )
{
  split_partition.clear();
  for( const auto &k:partition){
    split_blocks( k,  split_partition ) ;
    }

  assert( std::accumulate(split_partition.cbegin(), split_partition.cend(),0) <= std::accumulate(partition.cbegin(), partition.cend(),0));
}
\end{Piton}


@<separate entire partition@>=@#
void split_groups( const std::vector< unsigned int>& partition, std::vector<unsigned int>& split_partition )
{
@#
  split_partition.clear();
  @#
  for( const auto &k:partition){
  /* \newline \aths{\S~\ref{sec:split_blocks}} */
    split_blocks( k,  split_partition ) ;
    }
@#
  assert( std::accumulate(split_partition.cbegin(), split_partition.cend(),0) <= std::accumulate(partition.cbegin(), partition.cend(),0));
}



@*1 {\bf merge blocks}.
\label{sec:update_tree}


merge blocks and record continuing ones


\begin{Piton}
static void update_tree( std::vector<unsigned int>& tree, const std::vector<unsigned int>& mergers )
{
  assert(  mergers.size() > 0 );
  
  assert( static_cast<std::size_t>(std::accumulate(mergers.cbegin(), mergers.cend(), 0)) <= tree.size() ) ;
  std::shuffle( tree.begin(), tree.end(), rng) ;
  std::vector<unsigned int> newblocks (mergers.size());
  std::size_t j {} ;

  for( const auto &m: mergers ){
    newblocks[j] = std::accumulate( std::crbegin(tree), std::crbegin(tree) + m, 0);
    tree.resize( tree.size() - m ) ;
    ++j ;
  }
  tree.reserve( tree.size() + newblocks.size() );
  tree.insert( tree.end(), newblocks.cbegin(), newblocks.cend()  );
}
\end{Piton}


@<update tree@>=@#
static void update_tree( std::vector<unsigned int>& tree, const std::vector<unsigned int>& mergers )
{
@#
  assert(  mergers.size() > 0 );
@#  
  assert( static_cast<std::size_t>(std::accumulate(mergers.cbegin(),  mergers.cend(), 0)) <= tree.size() ) ;
@#

  std::shuffle( tree.begin(), tree.end(), rng) ;

@#

  std::vector<unsigned int> newblocks (mergers.size());

@#

  std::size_t j {} ;

@#

  for( const auto &m: mergers ){
    newblocks[j] = std::accumulate( std::crbegin(tree), std::crbegin(tree) + m, 0);
    tree.resize( tree.size() - m ) ;
    ++j ;
  }
  tree.reserve( tree.size() + newblocks.size() );
  tree.insert( tree.end(), newblocks.cbegin(), newblocks.cend()  );
}


@*1 {\bf sample time and partitions until a merger occurs}.
\label{sec:until_merger}


\begin{Piton}
static double until_merger(const std::size_t current_number_blocks, const std::vector<double>& v_lambdan, const std::vector<double>& v__cmf, std::vector<unsigned int>& v_merger_sizes, double &otimi)
{
  std::vector<unsigned int> v_partition  (SAMPLE_SIZE/2);
  v_partition.reserve(SAMPLE_SIZE/2) ;

  unsigned int lina {} ;
  double t {};

  auto newtime  = [](const double l, const double otime){
  return (RHO > DBL_EPSILON ? log1p( -(RHO * exp(- RHO * otime) / l) * log(gsl_rng_uniform_pos(rngtype)))/RHO : gsl_ran_exponential(rngtype, 1./l)) ;
  };

  auto samplemerger = [&v__cmf](void){
  
  unsigned int j {} ;
  const double u = gsl_rng_uniform( rngtype);
  
  while( u > v__cmf[j]){ ++j; } 
  
  return j ;
  } ;

  v_merger_sizes.clear() ;
  while( std::all_of(v_merger_sizes.cbegin(),  v_merger_sizes.cend(), [](const auto m){ return m < 2;}))
    {
      t += newtime( v_lambdan[current_number_blocks], otimi);
      otimi += t ;
      lina = samplemerger();
      readmergersizes( current_number_blocks, 1 + lina, v_partition) ;
      split_groups( v_partition, v_merger_sizes) ;
    }
  
  assert( v_merger_sizes.size() > 0) ;
  assert( std::all_of(v_merger_sizes.cbegin(),  v_merger_sizes.cend(), [](const auto m){ return m > 1;})); 
  return t ;
}
\end{Piton}


@<until merger@>=@#
static double until_merger(const std::size_t current_number_blocks, const std::vector<double>& v_lambdan, const std::vector<double>& v__cmf, std::vector<unsigned int>& v_merger_sizes, double &otimi)
{
@#
  std::vector<unsigned int> v_partition  (SAMPLE_SIZE/2);
  v_partition.reserve(SAMPLE_SIZE/2) ;
@#
  unsigned int lina {} ;
  double t {};

/* \newline \aths{new time $t =  \log\svigi{1 - \log(1-U)/\phi}/\rho$,
$\phi = \lambda_{n}\exp(\rho s)/\rho$, $s$ is cumulative time} */
  auto newtime  = [](const double l, const double otime){
  return (RHO > DBL_EPSILON ? log1p( -(RHO * exp(- RHO * otime) / l) * log(gsl_rng_uniform_pos(rngtype)))/RHO : gsl_ran_exponential(rngtype, 1./l)) ;
  };

@#

  auto samplemerger = [&v__cmf](void){
  @#
  unsigned int j {} ;
  const double u = gsl_rng_uniform( rngtype);
  @#
  while( u > v__cmf[j]){ ++j; } 
  @#
  return j ;
  } ;

@#

  v_merger_sizes.clear() ;
  while( std::all_of(v_merger_sizes.cbegin(),  v_merger_sizes.cend(), [](const auto m){ return m < 2;}))
    {
       t += newtime( v_lambdan[current_number_blocks], otimi);
      
      otimi += t ;
      @#
      lina = samplemerger();
         /* \newline \aths{\S~\ref{sec:readmergersizes}} */
       readmergersizes( current_number_blocks, 1 + lina, v_partition) ;
         /* \newline \aths{\S~\ref{sec:split_groups}} */
       split_groups( v_partition, v_merger_sizes) ;
    }
  @#
  assert( v_merger_sizes.size() > 0) ;
@#
  assert( std::all_of(v_merger_sizes.cbegin(),  v_merger_sizes.cend(), [](const auto m){ return m > 1;})); 

@#
return t ;
}


@*1 {\bf one experiment}.
\label{sec:one_experiment}


\begin{Piton}
static void one_experiment( std::vector<double>& v_l, std::vector<double>& v_r, const std::vector<double>& v__lambdan, const std::vector< std::vector< double > >&  a__cmf)
{

  std::vector<unsigned int> tree (SAMPLE_SIZE, 1) ;
  tree.reserve(SAMPLE_SIZE) ;
  std::fill(v_l.begin(), v_l.end(),0);
  
  double t {};
  double ot {} ;

  std::vector<unsigned int> v__merger_sizes {};
  v__merger_sizes.clear() ;
  v__merger_sizes.reserve(2*SAMPLE_SIZE) ;
  while( tree.size() > 1)
    {
      
      t = until_merger(tree.size(), v__lambdan, a__cmf[tree.size()], v__merger_sizes, ot);
      update_lengths(tree, v_l, t);
      update_tree(tree, v__merger_sizes);
    }
    @#
  assert( tree.back() == SAMPLE_SIZE);
  @#
  update_ri( v_r, v_l);
}
\end{Piton}


@<a single experiment@>=@#
static void one_experiment( std::vector<double>& v_l, std::vector<double>& v_r, const std::vector<double>& v__lambdan, const std::vector< std::vector< double > >&  a__cmf)
{
@#
  std::vector<unsigned int> tree (SAMPLE_SIZE, 1) ;
  tree.reserve(SAMPLE_SIZE) ;
  std::fill(v_l.begin(), v_l.end(),0);
  @#
  double t {};
  double ot {} ;
@#
  std::vector<unsigned int> v__merger_sizes {};
  v__merger_sizes.clear() ;
  v__merger_sizes.reserve(2*SAMPLE_SIZE) ;
@#

  auto update_lengths = [&tree, &v_l](const double t)
  {
  for (const auto &b: tree){
  v_l[0] += t;
  v_l[b] += t; }
  } ;

   while( tree.size() > 1)
    {
      /* \newline \aths{\S~\ref{sec:until_merger}} */     
      t = until_merger(tree.size(), v__lambdan, a__cmf[tree.size()], v__merger_sizes, ot);
      update_lengths( t);
      update_tree(tree, v__merger_sizes);
    }
    @#
  assert( tree.back() == SAMPLE_SIZE);
@#
   const double d = v_l[0];
   @#
   std::transform( v_l.cbegin(), v_l.cend(), v_r.begin(), v_r.begin(), [&d](const auto &x, const auto &y){return y + (x/d);});
}



@*1 {\bf approximate $\EE{R_{i}(n)}$}.
\label{sec:approximate_eri}


approximate $\EE{R_{i}(n)}$


\begin{Piton}
static void approximate_eri()
{
  std::vector<double> vri (SAMPLE_SIZE) ;
  vri.reserve(SAMPLE_SIZE) ;

  std::vector<double> v__l (SAMPLE_SIZE);
  v__l.reserve(SAMPLE_SIZE) ;

  std::vector<double> v__lambdan (SAMPLE_SIZE + 1) ;
  v__lambdan.reserve(SAMPLE_SIZE + 1) ;

  std::vector< std::vector< double > > a__cmfs (SAMPLE_SIZE + 1, std::vector<double> {} ) ;

  all_partitions_rates(v__lambdan, a__cmfs );

  int r = EXPERIMENTS + 1 ;

  while( --r > 0)
    {
      one_experiment(v__l,  vri, v__lambdan,  a__cmfs);
    }

  for( std::size_t i = 1; i < SAMPLE_SIZE; ++i){
    std::cout << (log( vri[i]/static_cast<double>(EXPERIMENTS)) - log1p(- vri[i]/static_cast<double>(EXPERIMENTS))) << '\n'; }
}
\end{Piton}


@<go ahead -- approximate $\EE{R_{i}(n)}$@>=@#
static void approximate_eri()
{
@#
  std::vector<double> vri (SAMPLE_SIZE) ;
  vri.reserve(SAMPLE_SIZE) ;

  std::vector<double> v__l (SAMPLE_SIZE);
  v__l.reserve(SAMPLE_SIZE) ;

  std::vector<double> v__lambdan (SAMPLE_SIZE + 1) ;
  v__lambdan.reserve(SAMPLE_SIZE + 1) ;

  std::vector< std::vector< double > > a__cmfs (SAMPLE_SIZE + 1, std::vector<double> {} ) ;

/* \newline \aths{\S~\ref{sec:all_partitions_rates}} */
  all_partitions_rates(v__lambdan, a__cmfs );

  int r = EXPERIMENTS + 1 ;

@#
  while( --r > 0)
    {/* \newline \aths{\S~\ref{sec:one_experiment}} */
      one_experiment(v__l,  vri, v__lambdan,  a__cmfs);
    }
@#
  for( std::size_t i = 1; i < SAMPLE_SIZE; ++i){
    std::cout << (log( vri[i]/static_cast<double>(EXPERIMENTS)) - log1p(- vri[i]/static_cast<double>(EXPERIMENTS))) << '\n'; }
}




@*1 {\bf main}. 
\label{sec:main}


the |main| module


@C


/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:veldi}} */
@<$e^{x}$@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<$\lambda_{n;k_{1},\ldots,k_{r};s}$ \eqref{eq:1}@>@#
/* \newline \aths{\S~\ref{sec:GenPartitions}} */
@<sizes of groups@>@#
/* \newline \aths{\S~\ref{sec:allmergers_sum_m}} */
@<sizes up to number@>@#
/* \newline \aths{\S~\ref{sec:ratesmergersfile}} */
@<partitions to files@>@#
/* \newline \aths{\S~\ref{sec:allmergers_when_n_blocks}} */
@<all partitions given number of blocks@>@#
/* \newline \aths{\S~\ref{sec:all_partitions_rates}} */
@<partitions generate all@>@#
/* \newline \aths{\S~\ref{sec:readmergersizes}} */
@<partitions read in@>@#
/* \newline \aths{\S~\ref{sec:split_blocks}} */
@<split partition into boxes@>@#
/* \newline \aths{\S~\ref{sec:split_groups}} */
@<separate entire partition@>@#
/* \newline \aths{\S~\ref{sec:update_tree}} */
@<update tree@>@#
/* \newline \aths{\S~\ref{sec:until_merger}} */
@<until merger@>@#
/* \newline \aths{\S~\ref{sec:one_experiment}} */
@<a single experiment@>@#
/* \newline \aths{\S~\ref{sec:approximate_eri}} */
@<go ahead -- approximate $\EE{R_{i}(n)}$@>@#

int main(int argc, const char * argv[])
{


  /* \newline \aths{\S~\ref{sec:rngs}} */
  setup_rng( static_cast<std::size_t>(atoi(argv[1])));
  /* \newline \aths{\S~\ref{sec:approximate_eri}} */
  approximate_eri();
  @#


  gsl_rng_free(rngtype) ;
@#

return 0 ;
}




@* {\bf  conclusions and bibliography}.
\label{sec:bib}

We approximate $\EE{R_{i}(n)}$ when the coalescent is the
time-changed
$\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent with
time-change function $G(t) = \int_{0}^{t}e^{\rho s}ds$ for a given
fixed $\rho \ge 0$.    
Figure~\ref{fig:graph} example approximation of 
$\EE{R_{i}(n)}$ 





\begin{SCfigure}[0.8][htb]
    \begin{tikzpicture}
      \begin{axis}[
        xlabel = $\log(i/n) - \log(1 - i/n)$,
        axis line style = {draw=none},
        tick style = {draw=none},
        xticklabels={draw=none},
        yticklabels={draw=none},
        legend pos=north east]
        \addplot table[col sep=comma] {forplottingfile1};
        %\addplot table[col sep=comma] {forplottingfile2};
        %\addplot table[col sep=comma] {forplottingfile3};
        %\addlegendentry{$\alpha = 0.05$}
        %\addlegendentry{$\alpha = 1$}
        %\addlegendentry{$\alpha = 1.5$}
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{ \it  An example approximation of $\EE{R_{i}(n)}$
 for the  given  parameter values and graphed as logits against
 $\log(i/n) - \log(1 - i/n)$ for $i = 1,2,\ldots, n-1$ where $n$ is
 sample size}}
      \label{fig:graph}
       \end{SCfigure}








\bibliographystyle{plain}


\begin{thebibliography}{99}

\bibitem[D2024]{D2024}
 Diamantidis,  Dimitrios and Fan,  Wai-Tong (Louis) and Birkner,
 Matthias and Wakeley,  John.  {Bursts of coalescence within
 population pedigrees whenever big families occur}.  Genetics Volume 227,  February
  2024.
 \\
 \url{https://dx.doi.org/10.1093/genetics/iyae030}.


\bibitem[CDEH25]{chetwyn-diggle_beta}
JA~Chetwyn-Diggle, Bjarki Eldon, and Matthias Hammer.
\newblock Beta-coalescents when sample size is large.
\newblock In preparation, 2025+.

\bibitem[DK99]{DK99}
P~Donnelly and T~G Kurtz.
\newblock Particle representations for measure-valued population models.
\newblock {\em Ann Probab}, 27:166--205, 1999.
\\
\url{https://dx.doi.org/10.1214/aop/1022677258}

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.


\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.



\bibitem[Pit99]{P99}
J~Pitman.
\newblock Coalescents with multiple collisions.
\newblock {\em Ann Probab}, 27:1870--1902, 1999.
\\
\url{https://dx.doi.org/10.1214/aop/1022874819}


\bibitem[Sag99]{S99}
S~Sagitov.
\newblock The general coalescent with asynchronous mergers of ancestral lines.
\newblock {\em J Appl Probab}, 36:1116--1125, 1999.
\\
\url{https://doi.org/10.1239/jap/1032374759}


\bibitem[Sch03]{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.
\\
\url{https://doi.org/10.1016/S0304-4149(03)00028-0}


\bibitem[S00]{S00}
J~Schweinsberg.
\newblock Coalescents with simultaneous multiple collisions.
\newblock {\em Electronic Journal of Probability}, 5:1--50, 2000.
\\
\url{https://dx.doi.org/10.1214/EJP.v5-68}



\bibitem[Tan11]{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}


@
\end{document}
