# This runs for each instance of terminal run by user leftynm
export PATH=$PATH:/usr/local/bin:/opt/local/bin

# This sets the gnuplot font default to something reliable
export GDFONTPATH=/Library/Fonts
export GNUPLOT_DEFAULT_GDFONT="Arial Narrow Bold.ttf"

# This controls the format of the terminal prompt
PS1="\h:\W \u:"

echo Welcome To The Machine
echo
cd Desktop

alias ..='cd ..'
alias ls='ls -F' 
alias ll='ls -lFh'
alias som='cd /Users/leftynm/Desktop/Research/ftran/SOM/Version3'
alias anis='cd /Users/leftynm/Desktop/Research/ftran/Anisotropic'
alias mix='cd /Users/leftynm/Desktop/Research/ftran/Mixed-Cells'
alias nyx='ssh leftynm@nyx-login.engin.umich.edu'
alias idl='/Applications/itt/idl/idl80/bin/idl'
alias blklght='ssh blacklight.psc.teragrid.org'
#alias octave='/opt/local/bin/octave'
#alias gfortran='/opt/local/bin/gfortran-mp-4.6'


echo Now activating the following aliases:
echo

echo ..='cd ..'
echo ls='ls -F'  
echo ll='ls -lFh'
echo idl='/Applications/itt/idl/idl80/bin/idl'
echo anis='cd /Users/leftynm/Desktop/Research/ftran/Anisotropic'
echo som='cd /Users/leftynm/Desktop/Research/ftran/SOM/Version3'
echo mix='cd /Users/leftynm/Desktop/Research/ftran/Mixed-Cells'
echo nyx='ssh leftynm@nyx-login.engin.umich.edu'
echo blklght='ssh blacklight.psc.teragrid.org'
#echo octave='/opt/local/bin/octave'
#echo gfortran='/opt/local/bin/gfortran-mp-4.6'

echo

# These commands make a pointer to gfortran 
#cd /usr/bin
#sudo ln -s /opt/local/bin/gfortran-mp-4.6 gfortran 
