import plot as plt
import generate_mesh as gm
import utilityFuncs as uf
import read_data as rd

def make_data( name, ndime, npoin, nelem, 
	       nnode, ngaus, nstre, lohi, 
	       coord, lnods, mass, v, stress, ep, force, stran, act, lelms, cell, keysg, keysl, keysp  ):
	dataFile = name
	stObj = st.makeDATA( write = open( dataFile, 'w' ) )
	d = dict( ndime = ndime, npoin = npoin, nelem = nelem, 
		  nnode = nnode, ngaus = ngaus, nstre = nstre, 
		  lohi  = lohi , coord = coord, lnods = lnods, mass =  mass, vel = v, sigma = stress, ep = ep, force = force, strain = stran, act = act, lelms = lelms, mn = cell, keysp = keysp, keysg = keysg, keysl = keysl )
	stObj.femPyData( **d )
	print 'data.txt created!'

def SampleDistribution( epsilon, gmodu, epsilonRecovery,
			meanGmodu, meanYieldStrain, meanRecoveryStrain, width, crit, FRICT, mnunx, resid, gauss = None, uniform = True ):
	nelem = len( gmodu )
	# --- draw from uniform distribution
#	meanE = 0.5 * meanGmodu * meanYieldStrain * meanYieldStrain
#	lamda = 700.0 
#	minE  = 0.25 * meanGmodu * 0.035 * 0.035 #meanRecoveryStrain * meanRecoveryStrain 
#	minE  = 0.5 * meanGmodu * 0.035 * 0.035 #meanRecoveryStrain * meanRecoveryStrain 
#	meane = 0.5 * meanGmodu * ( 0.1 ) ** 2
#	if crit == 'mohr-col':
#		meane = 0.5 * meanGmodu * ( 0.1 / cos( FRICT ) ) ** 2
	lamda = 1.0 / ( mnunx - resid )
	print 'lamda = %e' %lamda
	print 'minE = %e'%resid
	yieldStrain = []
	for ielem in xrange( nelem ):
		x = rnd.expovariate( lamda ) 
#		x = rnd.random() * 2 * ( mnunx - resid ) + resid
#		yieldStrain.append( 2.0 * ( x + minE ) ** 0.5 )
#		yieldStrain.append( ( 2.0 * ( x + minE ) ) ** 0.5 )
		yieldStrain.append( x + resid )
#		yieldStrain.append( x )
#		x = rnd.random()
#		yieldStrain.append( 1.0e+08 )
	epsilon[ : ] = yieldStrain[ : ] # --- update epsilon
	epsilonRecovery[ : ] = [ meanRecoveryStrain for i in yieldStrain ] # --- update epsilonRecovery:
	# --- shear moduli
	gmodu[ : ] = [ 1.0 for i in xrange( nelem ) ] # --- update gmodu

#	sfile=open('/scratch/kkarimi/myRuns/20x20/test/Run0/mohr/uniax.xyz')
#	natom=9646
#	[sfile.readline() for i in xrange(9)]
#	epsilon[:]=[float(sfile.readline().split()[1]) for i in xrange(natom)]

def get_mod( name, lnods, K, Gbar, ystrn, recovery_strain, FRICT, crit, resid, reduc, mnunx, pcrit, rsidf, rsidp, redup, reduf ):
# --- create moduli.txt
	modFile = name
	nelem = len( lnods )
	File = open( modFile, 'w' )
	print >> File, '#type, kbulk, gmod, uniax, epsrc, frict, hards, min uniax, min frict, reduc, redup, reduf, <uniax>, pcrit, min pcrit'

	epsilon = [ 0.0 for ielem in xrange( nelem ) ]
	gmodu = [ 0.0 for ielem in xrange( nelem ) ]
	epsilonRecovery = [ 0.0 for ielem in xrange( nelem ) ]


	SampleDistribution( epsilon, gmodu, epsilonRecovery , # --- output
			    Gbar, ystrn, recovery_strain, width = 0.25, gauss = True, uniform = None, crit = crit, FRICT = FRICT, mnunx = mnunx, resid = resid ) # --- width
	meanFricAngle = ( FRICT )  #--- mean friction angle in mohr-coulomb criterion
#	meanFricAngle = 80.0 * pi / 180.0 #--- mean friction angle in mohr-coulomb criterion
#	fricAngle = []
#	for ielem in xrange( nelem ):
#		fricAngle.append( rnd.expovariate( 1.0 / meanFricAngle ) )
	fricAngle = nelem * [ meanFricAngle ] 
	'''
	# --- draw histogram
	f = open( 'pdf.txt','w')
	hmObj = hm.histogram( 25, min( gmodu ), max( gmodu ) ) # --- gmodu
	for item in gmodu:
		hmObj.whichBin( item, [ 0.0 ] )
	for i in hmObj.res():
        	for j in i: print >> f, j, '\t',
        	print >> f
        print >> f
	hmObj = hm.histogram( 25, min( epsilon ), max( epsilon ) ) # --- yield strain
	for item in epsilon:
		hmObj.whichBin( item, [ 0.0 ] )
	for i in hmObj.res():
        	for j in i: print >> f, j, '\t',
        	print >> f
	f.close()
	'''
	# --- output
#	assert min( epsilon ) > maxPlaStran

	for ielem in xrange( nelem ):
#		epsilon[ ielem ] = 1.0e+03
#		if ielem != 64:
#			epsilon[ ielem ] = 1.0e+03
#			epsilon[ ielem ] = 1.5 * epsilon[ 64 ] * gmodu[ 64 ] #--- shear band
		print >> File, ielem, gmodu[ ielem ] * K / Gbar, gmodu[ ielem ], epsilon[ ielem ] * gmodu[ ielem ], epsilonRecovery[ ielem ], fricAngle[ ielem ], hards, resid, rsidf, reduc, redup, reduf, mnunx, pcrit, rsidp
#		if i == 0:
#			print >> File, i, K, Gbar, ystrn * Gbar * 0.1, recovery_strain
#		else:
#			print >> File, i, K, Gbar, ystrn * Gbar, recovery_strain
	File.close()
	print '%s created!' %modFile
# -----------------------------------------------------------------------------------


def get_lnods_triangles( alpha, lohi, coord, nnode ):
	gmObj = gm.mesh2D( lohi = lohi,
		      coord = coord, nnode = nnode )
	lnods, elcod = gmObj.get_lnods( alpha )
	return	lnods, elcod

def get_coord_lnods( nnode, rcut, file ):

	# --- triangular elements
	if 1:
# --- read lammps output 
		someFile = open( file )
		for i in xrange( 2 ):
			someFile.readline()
		npoin = int(someFile.readline().split()[0])
		for i in xrange( 11 ):
			someFile.readline()
		[ xlo, xhi ] = map( float, someFile.readline().split()[0:2] )
		[ ylo, yhi ] = map( float, someFile.readline().split()[0:2] )
#		print [ ylo, yhi ]
		for i in xrange( 7 ):
			someFile.readline()
		tmp_list = [ map( float, someFile.readline().split()[ 2 : 4 ] ) for i in xrange(npoin) ]
		lx = xhi - xlo
		ly = yhi - ylo
#		rcut = [xhi - xlo,yhi-ylo] #10.0 # --- cut a square of length rcut 
		xmin = xlo #0.5 * ( xlo + xhi - rcut[ 0 ] )
		ymin = ylo #0.5 * ( ylo + yhi - rcut[ 1 ] )
		xmax = xhi #xmin + rcut[ 0 ]
		ymax = yhi #ymin + rcut[ 1 ]
		tmp = []
		for item in tmp_list: #coord:
			[ x, y ] = item
#			assert xmin <= x < xmax, '%s <= %s < %s' %(xmin, x, xmax) 
			if not xmin <= x < xmax:
				x = ( x - xmin ) % lx + xmin
				assert xmin <= x < xmax, '%s <= %s < %s' %(xmin, x, xmax) 
#			assert ymin <= y < ymax, '%s <= %s < %s' %(ymin, y, ymax) 
			if not ymin <= y < ymax:
				y = ( y - ymin ) % ly + ymin
			tmp.append( [x,y] )
#			if xmin <= x < xmax and ymin <= y < ymax: 
#				tmp.append( [ ( item[ 0 ] - xmin ), ( item[ 1 ] - ymin ) ] )
		coord = {}
		for index, items in zip( xrange( sys.maxint ), tmp ):
			coord[ index ] = items[ : ]
#		coord = [ [ ( i[ 0 ] - xlo ) / lx - 0.5, ( i[ 1 ] - ylo ) / ly - 0.5 ] for i in coord ] # --- rescale xyz
		npoin = len( coord )
#		lx = ly = 1.0
#		coord = [ [ 0.0, 0.0 ] for i in xrange( m ) for j in xrange( n ) ]
#		for i in xrange( - ( m / 2 ), m / 2 + m % 2 ):
#			y = i * ly / m + ( rnd.random() - 0.5 ) * 0.25 / n 
#			for j in xrange( - ( n / 2 ), n / 2 + n % 2 ):
#				x = j * lx / n + 0.5 * y * xy / ly + ( rnd.random() - 0.5 ) * 0.25 / n	
#				index = ( i % m ) * n + j % n
#				coord[ index ] = [ x, y ]
#		coord = [ [ rnd.random() - 0.5 for idime in xrange( ndime ) ] for ipoin in xrange( npoin ) ]
#		xlo = 0.0
#		xhi = rcut[ 0 ] #1.0
#		ylo = 0.0
#		yhi = rcut[ 1 ]
		lohi = [ [ xlo, xhi ], [ ylo, yhi ], 0.0 ]
		lx = xhi - xlo
		ly = yhi - ylo

		lnods, elcod = get_lnods_triangles( alpha = 0.6, lohi = lohi, coord = coord, nnode = 3 )

	return	coord, lnods, elcod, lohi

def plot_mesh( lnods, elcod, nnode, coord, lohi, name, text = None ):

	# --- plot
	pltObj = plt.plot( coord, lohi, name )
	pltObj.mesh( elcod )
	[ [ xlo, xhi ], [ ylo, yhi ], xy ] = lohi
	lx = xhi - xlo
	ly = yhi - ylo
#	h = [ [ lx, xy ], [ 0.0, ly ] ]
	i = 0
	j = 0
	N=1
	M=1
	h = [ [ ( lx / N ), 0.0 ], [ 0.0, ( ly / M ) ] ]
	center = mat.matrix( h ) * [ 0.5, 0.5 ]
#	center[ 0 ] += xlo
#	center[ 1 ] += ylo
	center[ 0 ] += j * ( lx / N )
	center[ 1 ] += i * ( ly / M )
	pltObj.plotBox( h = h, center = center )
	if text:
		pltObj.writeText( xy_list = text[ 0 ], str_list = text[ 1 ], size = 0.5 )
	pltObj.get()
	print 'open %s.pdf' %name
	#---
class junk( rd.read_data ):
	def __init__( self ):
		self.npoin = len( coord )
		self.ndime = len( coord[ 0 ] )
		self.nnode = 3 # --- number of nodes
		self.nsvab = self.npoin * self.ndime
		self.nevab = self.nnode * self.ndime
		self.lohi = lohi
		self.coord = coord
		self.ngaus = 1 # --- gauss points
		self.nstre = 3 # length of stress vector
		self.nelem = len( lnods )
		self.mass = [ 0.0 for ielem in xrange( self.nelem ) ]
		for ielem in xrange( self.nelem ):
			self.mass[ ielem ] = self.npoin * 1.0 #rho[ iRun ] # --- density
		self.lnods = lnods
# -----------------------------------------------------------------------------------

if __name__ == '__main__':
	import random as rnd
	import scipy.spatial as sp
	from math import *
	import sys
	import os
	import numpy as np

	import expand_periodic_box as pbc
	import matrix as mat
	import dist_ij as dij
	import rmDuplicates as rm
	import mapping as mp
	import stream as st
	import compute_mass_matrix as cm
	import histogram as hm

#--- 
	pathForDump = sys.argv[ 1 ] #--- path for a discrete system (dump.xyz )
	boxDim = [ [ float( sys.argv[ 2 ] ), float( sys.argv[ 3 ] ) ] ]
	maxPlaStran = float( sys.argv[ 3 ] ) #0.01 # --- read max plastic strain
	hards = float( sys.argv[ 4 ] ) #- 1.9 # --- H=-0.2\mu shouldn't be -1.0 for von-mises plasticity
	fyild = float( sys.argv[ 5 ] ) #-1.0e-04 #--- distance to sy
	K = float( sys.argv[ 6 ] ) #2.0
	Gbar = float( sys.argv[ 7 ] ) #1.0
	pcrit = K * float( sys.argv[ 8 ] ) #1.0e+06 #--- press cap only if press-dependent (should be big for no dilatancy)
	rsidp = float( sys.argv[ 9 ] ) * pcrit #0.0 * pcrit #--- residual press (weakening)
	redup = float( sys.argv[ 10 ] ) #0.0 #--- reduction in cap press
	FRICT = float( sys.argv[ 11 ] ) * pi / 180 # 65.0 * pi / 180 #--- deg
	reduf = float( sys.argv[ 12 ] ) #0.0 #--- reduction in friction
#	Reduc = [ 0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.99 ]
#	RSIDF = [ ( 1.0 - i ) * FRICT for i in Reduc ]
#	frict = [ 35.0, 45.0, 55.0, 65.0, 75.0 ]
#	frict = [ items * pi / 180.0 for items in frict ]
	rsidf = float( sys.argv[ 13 ] ) * FRICT  #1.0 * FRICT #--- residual friction (weakening)
	mnunx = float( sys.argv[ 14 ] ) #0.1 #--- mean yield stress
	reduc = float( sys.argv[ 15 ] ) #0.0 #1.0 #--- reduction in yield stress for damage model
	resid = float( sys.argv[ 16 ] ) * mnunx #0.1 * mnunx #0.0 * mnunx #--- residual yield stress (weakening)
#	Reduc = [ 0.0, 0.1, 0.2, 0.4, 0.8 ]
#---
	pinit = float( sys.argv[ 17 ] ) * mnunx #8.0 * mnunx #1.0 * pcrit #--- initial press 
#	PINIT = [ 2.0, 4.0, 8.0, 16 ]
#	PINIT = [ items * mnunx for items in PINIT ]
#---
	UNFRM = bool( int( sys.argv[ 18 ] ) ) #True #--- uniform initial stress
	ZSHER = bool( int( sys.argv[ 19 ] ) ) #None #--- zero shear stress
#---
	crit = sys.argv[ 20 ] #'mohr-col' 'von-mises'
	shear = sys.argv[ 21 ]# 'pure', 'simple' 
	rcut = float( sys.argv[ 22 ] ) #4.0 # --- cell size
#---
	if 1:
		index = 0
		for iRun in xrange( 1 ): #nRun ):
#			reduf = Reduc[ iRun ]
#			rsidf = RSIDF[ iRun ]
			coord, lnods, elcod, lohi = get_coord_lnods( nnode = 3, rcut = boxDim[ 0 ], file = pathForDump ) # --- generates the mesh
#			plot_mesh( lnods = lnods, nnode = 3, 
#				       coord = coord, lohi = lohi, elcod = elcod, name = 'mesh' ) #, text = someList )
			[ [ xlo, xhi ], [ ylo, yhi ], xy ] = lohi
			print 'lohi=%s'%lohi
			assert min( [ coord[ xy ][ 0 ] for xy in coord ] ) >= xlo
			assert max( [ coord[ xy ][ 0 ] for xy in coord ] ) <= xhi
			assert min( [ coord[ xy ][ 1 ] for xy in coord ] ) >= ylo
			assert max( [ coord[ xy ][ 1 ] for xy in coord ] ) <= yhi
			get_mod( K = K, Gbar = Gbar, ystrn = 0.1, recovery_strain = maxPlaStran, lnods = lnods, name = 'moduli.txt', \
					FRICT = FRICT, crit = crit, resid = resid, rsidf = rsidf, rsidp = rsidp, \
					reduc = reduc, mnunx = mnunx, pcrit = pcrit, redup = redup, reduf = reduf ) # ---assign shear & bulk moduli: check func get_mod & moduli.h 
			print 'mean = uniax', np.mean( [ float( line.split()[ 3 ] ) for line in open( 'moduli.txt' )  if line[ 0 ] != '#' ] )
			print 'min uniax = ', min( [ float( line.split()[ 3 ] ) for line in open( 'moduli.txt' )  if line[ 0 ] != '#' ] )
			print 'max uniax = ', max( [ float( line.split()[ 3 ] ) for line in open( 'moduli.txt' )  if line[ 0 ] != '#' ] )
			print 'mean angle = ', np.mean( [ float( line.split()[ 5 ] ) for line in open( 'moduli.txt' )  if line[ 0 ] != '#' ] ) * 180.0 / pi
		# --- compute mass matrix
			femObj = junk()
			cmObj = cm.compute_mass_matrix()
			cmObj.get_mass_matrix( femObj = femObj, diagonal_mass = True )
			mass = [ 0.0 for ipoin in xrange( femObj.npoin ) ]
			isvab = 0
			for ipoin in xrange( femObj.npoin ):
				for idime in xrange( femObj.ndime ):
					mass[ ipoin ] = femObj.tmass[ isvab, isvab ] / cmObj.vol #---mass
					isvab += 1;
			# --- initial vel	
			v = [ [ 0.0, 0.0 ] for i in xrange( len( coord ) ) ]
		# --- initial force	
			force = [ [ 0.0, 0.0 ] for i in xrange( len( coord ) ) ]
		# --- initial stress: add some noise to decorrelate elements 	
			uniax = min( [ float( line.split()[ 3 ] ) for line in open( 'moduli.txt' )  if line[ 0 ] != '#' ] ) #--- minimum yield stress
#			if crit == 'mohr-col':
#				uniax = min( [ pinit * sin( float( line.split()[ 5 ] ) ) + float( line.split()[ 3 ] ) * cos( float( line.split()[ 5 ] ) ) for line in open( 'moduli.%s.%s.txt' % ( index, iRun ) )  if line[ 0 ] != '#' ] )
#			print 'uniax = ', uniax
#			if ZSHER:
#				fyild = uniax
#			else:
#				assert uniax - fyild >= 0.0, 'reduce fyild!'
			preys = uniax
			if crit == 'mohr-col':
				preys = pinit * sin( FRICT ) + uniax * cos( FRICT )
			if shear == 'simple':
				if UNFRM:
					stress = [ [ 0.0, 0.0, 1.0 * ( preys + fyild ) ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
				else:			
					stress = [ [ 0.0, 0.0, rnd.random() * ( preys + fyild ) ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
				if crit == 'mohr-col':
					if UNFRM:
						stress = [ [ - pinit, - pinit, 1.0 * ( preys + fyild ) ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
					else:
						stress = [ [ - pinit, - pinit, rnd.random() * ( preys + fyild ) ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
#
			if shear == 'pure':
				if UNFRM:
					stress = [ [ 1.0 * ( preys + fyild ), - 1.0 * ( preys + fyild ), 0.0 ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
				else:
					Delta = [ rnd.random() * ( preys + fyild ) for i in xrange( len( lnods ) * femObj.ngaus ) ]
					stress = [ [ delta, - delta, 0.0 ] for delta in Delta ]
				if crit == 'mohr-col':
					if UNFRM:
						stress = [ [ - pinit + 1.0 * ( preys + fyild ), - pinit - 1.0 * ( preys + fyild ), 0.0 ] for i in xrange( len( lnods ) * femObj.ngaus ) ]					
					else:
#						Delta = [ ( preys + fyild ) for i in xrange( len( lnods ) * femObj.ngaus ) ]
						Delta = [ rnd.random() * ( preys + fyild ) for i in xrange( len( lnods ) * femObj.ngaus ) ]
#						mu = 0.75
#						Delta = [ rnd.gauss( mu, 0.1 * mu ) * ( preys + fyild ) for i in xrange( len( lnods ) * femObj.ngaus ) ]
						stress = [ [ - pinit + delta, - pinit - delta, 0.0 ] for delta in Delta ]
			if crit == 'mohr-col':
				print 'fymax=', max( [ 0.5*abs(stress[i-1][0]-stress[i-1][1])-( pinit * sin( float( line.split()[ 5 ] ) ) + float( line.split()[ 3 ] ) * cos( float( line.split()[ 5 ] ) )) for line, i in zip(open( 'moduli.txt' ),xrange(sys.maxint)) if line[ 0 ] != '#' ] )
			else:
				print 'fymax=', max( [ 0.5*(stress[i-1][0]-stress[i-1][1])- float( line.split()[ 3 ] ) for line, i in zip(open( 'moduli.txt' ),xrange(sys.maxint)) if line[ 0 ] != '#' ] )

			#---
			stran = [ [ 0.0, 0.0, 0.0 ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
			act = [ 0 for i in xrange( len( lnods ) * femObj.ngaus ) ]
	# --- initial plastic strain
			recovery_strain = maxPlaStran
			ep = [ [ 0.0, 0.0, 0.0 ] for i in xrange( len( lnods ) * femObj.ngaus ) ]
	# --- msd analysis: construct square cells and find triangular elements inside each cell
			lelms = {}
			N = int( ( xhi - xlo ) / rcut ) # --- number of cells in x
			M = int( ( yhi - ylo ) / rcut )
			nelem = len( lnods )
			for ielem in xrange( nelem ):
				for nodei in lnods[ ielem ]:
					xIndx = int( N * ( coord[ nodei ][ 0 ] - xlo ) / ( xhi - xlo ) ) % N  # --- periodic cell
					yIndx = int( M * ( coord[ nodei ][ 1 ] - ylo ) / ( yhi - ylo ) ) % M
					assert 0 <= xIndx < N
					assert 0 <= yIndx < M
					lelms.setdefault( yIndx * N + xIndx, [] ).append( ielem ) # --- key: cell index, val = element id
			# --- rm repetative elements
			tmp = {}
			for icell in xrange( M * N ):
				e = {}
				for item in lelms[ icell ]:
					e.setdefault( item, [] ).append( 0 )
				someList = e.keys()
				someList.sort()
				assert len( someList ) > 1
				tmp[ icell ] = someList
			lelms = tmp
#			print 'elements', len( lelms[ 0 * N + 0 ] ), lelms[ 0 * N + 0 ]
			'''
			# --- plot 
			positions = []
			string_list = []
			for item, Index in zip( elcod, xrange( sys.maxint) ):
				xc = 0.0
				yc = 0.0
				for i in item:
					xc += i[ 0 ]
					yc += i[ 1 ]
				positions.append( [ xc / 3.0, yc / 3.0 ] )
			 	string_list.append( str( Index ) )
			text = [ positions, string_list ]
			plot_mesh( lnods = lnods, nnode = 3, 
				    coord = coord, lohi = lohi, elcod = elcod, name = 'mesh%s' % index, text = text )
			'''
	# --- output data.txt
			keysp = [ i for i in xrange( len( coord ) )]
			keysg = [ i for i in xrange(femObj.ngaus*len( lnods ))]
			keysl = [ i for i in xrange(len( lnods ))]	
			make_data( name = 'data.txt', ndime = 2, npoin = len( coord ), nelem = len( lnods ), 
			       	   nnode = femObj.nnode, ngaus = femObj.ngaus, nstre = femObj.nstre, lohi = lohi, 
			           coord = coord, lnods = lnods, mass = mass, v = v, stress = stress, ep = ep, force = force, stran = stran, act = act, lelms = lelms, cell = [ M, N ], keysp = keysp, keysg = keysg, keysl = keysl  )

