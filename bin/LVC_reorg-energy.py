#!/usr/bin/env python2

from optparse import OptionParser

au2rcm = 219474.631

# =========================================================
def find_lines(nlines,match,strings):
  smatch=match.lower().split()
  iline=-1
  while True:
    iline+=1
    if iline==len(strings):
      return []
    line=strings[iline].lower().split()
    if tuple(line)==tuple(smatch):
      return strings[iline+1:iline+nlines+1]

# =========================================================

def read_V0_simple(SH2LVC, fname='V0.txt'):
  """"
  Reads information about the ground-state potential from V0.txt.
  Returns the displacement vector.
  """
  try:
    f=open(fname)
  except IOError:
    print 'Input file %s not found.'%fname
    sys.exit(20)
  v0=f.readlines()
  f.close()

  # Frequencies (a.u.)
  tmp = find_lines(1, 'Frequencies',v0)
  if tmp==[]:
    print 'No Frequencies defined in %s!'%fname
    sys.exit(22)
  SH2LVC['Om'] = [float(o) for o in tmp[0].split()]

  # Normal modes in mass-weighted coordinates
  tmp = find_lines(len(SH2LVC['Om']), 'Mass-weighted normal modes', v0)
  if tmp==[]:
    print 'No normal modes given in %s!'%fname
    sys.exit(23)
  SH2LVC['V']  = [map(float,line.split()) for line in tmp] # transformation matrix

# =========================================================

def read_SH2LVC_simple(fname='LVC.template'):
  # reads LVC.template, deletes comments and blank lines
  SH2LVC={}
  try:
    f=open(fname)
  except IOError:
    try:
      f=open('SH2LVC.inp')
    except IOError:
      print 'Input file "LVC.template" not found.'
      sys.exit(24)
  sh2lvc=f.readlines()
  f.close()

  disp = read_V0_simple(SH2LVC, sh2lvc[0].strip())

   # Add the vertical energies (epsilon)
  # Enter in separate lines as:
  # <n_epsilon>
  # <mult> <state> <epsilon>
  # <mult> <state> <epsilon>

  tmp = find_lines(1, 'epsilon',sh2lvc)
  if not tmp==[]:
    SH2LVC['eps'] = []
    neps = int(tmp[0])
    tmp = find_lines(neps+1, 'epsilon', sh2lvc)
    for line in tmp[1:]:
      words = line.split()
      SH2LVC['eps'].append((int(words[0])-1, int(words[1])-1, float(words[-1])))

  #for imult in range(nmult): print numpy.array(HMCH[imult])

  # Add the intrastate LVC constants (kappa)
  # Enter in separate lines as:
  # <n_kappa>
  # <mult> <state> <mode> <kappa>
  # <mult> <state> <mode> <kappa>

  tmp = find_lines(1, 'kappa', sh2lvc)
  if not tmp==[]:
    SH2LVC['kappa'] = []
    nkappa = int(tmp[0])
    tmp = find_lines(nkappa+1, 'kappa', sh2lvc)
    for line in tmp[1:]:
      words = line.split()
      SH2LVC['kappa'].append((int(words[0])-1, int(words[1])-1, int(words[2])-1, float(words[-1])))


  # Add the interstate LVC constants (lambda)
  # Enter in separate lines as:
  # <n_lambda>
  # <mult> <state1> <state2> <mode> <lambda>
  # <mult> <state1> <state2> <mode> <lambda>

  tmp = find_lines(1, 'lambda', sh2lvc)
  if not tmp==[]:
    SH2LVC['lam'] = []
    nlam = int(tmp[0])
    tmp = find_lines(nlam+1, 'lambda', sh2lvc)
    for line in tmp[1:]:
      words = line.split()
      SH2LVC['lam'].append((int(words[0])-1, int(words[1])-1, int(words[2])-1, int(words[3])-1, float(words[-1])))

  return SH2LVC

# ============================================================================

def setup_parser():
    usage='''
    %s [options]

    This script computes Huang-Rhys factors S_i and partial reorganisation
    energies lam_i from LVC parameters.
    
    S_i   = 0.5 * kappa_i^2 / omega_i^2
    lam_i = 0.5 * kappa_i^2 / omega_i'''%__file__.split('/')[-1]

    parser = OptionParser(usage=usage)
    parser.add_option('-m', dest='m', type=int, nargs=1, default=1, help="Multiplicity of target state (default=1)")
    parser.add_option('-s', dest='s', type=int, nargs=1, default=2, help="Index of target state (default=2)")
    parser.add_option('-p', dest='p', default=False, help="Plot using pylab", action="store_true")
    
    return parser
# ============================================================================

def plot_lam(Om, lam, pname='reorg.png'):
    import pylab
    from matplotlib import rc

    pylab.figure(figsize=(8,6))
    rc('font', size=14)

    for i, omval in enumerate(Om):
        orcm = omval * au2rcm
        pylab.plot([orcm, orcm], [-1., lam[i] * au2rcm], 'ro-')

    pylab.xlabel('Wavenumber (1/cm)')
    pylab.ylabel('Reorganisation en. (1/cm)')
    pylab.axis(ymin=-10.)
    pylab.savefig(pname)

    print("Plot file %s created."%pname)
# ============================================================================
def plot_S(Om, S, pname='Huang_Rhys.png'):
    import pylab
    from matplotlib import rc

    pylab.figure(figsize=(8,6))
    rc('font', size=14)

    for i, omval in enumerate(Om):
        orcm = omval * au2rcm
        pylab.plot([orcm, orcm], [-1., S[i]], 'ro-')

    pylab.xlabel('Wavenumber (1/cm)')
    pylab.ylabel('Huang-Rhys factor')
    pylab.axis(ymin=-.01)
    pylab.savefig(pname)

    print("Plot file %s created."%pname)
# ============================================================================

def main():
    (options, args) = setup_parser().parse_args()
    SH2LVC=read_SH2LVC_simple()

    jmult  = options.m - 1
    jstate = options.s - 1

    S   = [0. for omval in SH2LVC['Om']]
    lam = [0. for omval in SH2LVC['Om']]

    ncontrib = 0
    for k in SH2LVC['kappa']:
        (imult, istate, i, val) = k
        if imult==jmult and istate==jstate:
            ncontrib += 1
            S[i] =   0.5 * val**2 / SH2LVC['Om'][i]**2
            lam[i] = 0.5 * val**2 / SH2LVC['Om'][i]

    print "Total reorganisation energy: %10.2f cm-1"%(sum(lam) * au2rcm)
    print "Number of contributing modes: %i/%i\n"%(ncontrib, len(SH2LVC['Om']))

    print "%5s %10s %10s %10s"%('', 'om/cm-1', 'S_i', 'lam_i/cm-1')
    print 32*"-"
    for i, omval in enumerate(SH2LVC['Om']):
        print "%5i %10.2f %10.3f %10.2f"%(i+1, omval * au2rcm, S[i], lam[i] * au2rcm)

    print 32*"-"
    print "%16s %10.3f %10.2f"%('', sum(S), sum(lam) * au2rcm) 

    if options.p:
        plot_lam(SH2LVC['Om'], lam)
        plot_S(SH2LVC['Om'], S)

# ============================================================================
if __name__ == '__main__':
  main()
