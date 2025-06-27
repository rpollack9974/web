from sage.geometry.newton_polygon import NewtonPolygon

class ghost(SageObject):
	"""This class represents the ghost series for a prime p and a fixed tame level Gamma_0(N).
		By ghost series, we mean the (p-1)/2 power series for each even component of weight space

		To initialize a ghost series simply type:

			ghost(p,N) 

		or to specify the number of coefficients computed:

			ghost(p,N,num_coefs=??)
		
		The key functions of this class are:

		1) slopes
		2) wadic_slopes

		If g is the ghost series, the first is called by:

			g.slopes(k) or g.slopes(k,num=??) 

		and returns the ghost slopes in weight k.

		The second is called by:

			g.wadic_slopes(comp) or g.wadic_slopes(comp,num=??)

		and returns the w-adic slopes on the component comp.
		(Components are represented by even integers between 0 and p-2.)

		A fancier commend is:

		3) global_halo

		This draws a picture of the ghost spectral curve.  To call this function type:

			g.global_halo(comp,sheets,max_v)

		where comp is again the component, sheets is the number of 'sheets' in the picture, 
		and max_v is the maximum valuation of weights used.  For example, for p=2 and N=1,

			g.global_halo(0,10,10)

		will draw the first 10 sheets of the global halo over the region of weight space 
		for weights with (non-integral) valuation between 0 and 10.  

	 	There is also an optional argument 'central_wt' where one can change the weight we
	 	are expanding around.

		Lastly, to access the actual ghost zeroes type:

			g[comp][i]

		to see a list with multiplicities of the ghost series of the i-th coefficient
		on component comp.  What you see should be self-explanatory except if p=2 and N>1.
		Then the "modified" ghost zeroes are represented by negative integers.
	"""

	def compute_ghost_series(self,num_coefs=10):
		"""This function computes and stores the relevant ghost seris in 'shell' form --- that is,
			for each component comp, it sets self.series[comp] to a list whose i-th term is a 
			list of the ghost zeroes of the i-th coefficient with multiplicities.

			Inputs:
				num_coefs (optional) -- an integer representing how many coefficients are computed.
		"""
		p=self.p
		N=self.N

		ghost_coefs = [[[] for i in range(num_coefs+1)] for a in range(p-1)]

		## Starting at weight 2, we run through weights in the component,
		## compute the associated indices, and then record the weights at those
		## indices with appropriate multiplicities

		k=2;
		if k==0:
			k=k+p-1
		f1 = dimension_cusp_forms(N,k)
		f2 = dimension_cusp_forms(N*p,k)-2*f1

		inds = range(f1+1,f1+f2)
		while (len(inds)==0 or inds[0]<=num_coefs+1):
			## This loops adds the weights to the appropriate indices with the appropriate multiplicities
			for m in range((len(inds)+1)/2):
				if m < floor(len(inds)/2):
					if inds[m]<=num_coefs:
						ghost_coefs[k%(p-1)][inds[m]] += [(k,m+1)]
					if (inds[len(inds)-1-m]<=num_coefs):
						ghost_coefs[k%(p-1)][inds[len(inds)-1-m]] += [(k,m+1)]
				else:
					if inds[m]<=num_coefs:
						ghost_coefs[k%(p-1)][inds[m]] += [(k,m+1)]
			k = k+1
			f1 = dimension_cusp_forms(N,k)
			f2 = dimension_cusp_forms(N*p,k)-2*f1
			inds = range(f1+1,f1+f2)
		self.series=ghost_coefs
		if p==2 and N>1:
			self.form_p2_modification()

	def __init__(self,p,N,num_coefs=10):
		"""Initializes the ghost series

		Inputs:
			p - prime
			N - tame level
			num_coefs (optional) - number of coefficients computed
		"""
		self.p=p
		self.N=N
		self.compute_ghost_series(num_coefs=num_coefs)

	def __repr__(self):
		"""Returns the string <Ghost series for p and N>"""
		return "Ghost series for p="+str(self.p)+" and N="+str(self.N)

	def __getitem__(self,a):
		"""Returns the power series (in shell form) on a-th component"""
		return self.series[a]

	def ghost_zeroes(self,comp,i):
		"""Returns a list of the ghost series of the i-th coefficient on component comp"""
		return self.series[comp][i]

	def num_coefs(self,comp):
		"""Returns the number of coeffcients computed on component comp"""
		return len(self[comp])

	def slopes(self,k,num=None):
		"""Returns the slopes in weight k --- num-many or all computed if term==None"""

		NP = []
		p=self.p
		comp=k%(p-1)
		if num==None:
			d = self.num_coefs(comp)
		else:
			### HACKING HERE (UNCLEAR HOW MANY TERMS NEED TO BE USED TO GET CORRECT SLOPES)
			if self.num_coefs(comp)<num+10:
				self.compute_ghost_series(num_coefs=num+10)
			d = min(self.num_coefs(comp),num+10) 
		if p == 2:
			e = 2
		else:
			e = 1
		for i in range(d):
			y = 0
			for ss_wt in self[comp][i]:
				mult = ss_wt[1]
				k_ss = ss_wt[0]
				if k_ss == "p":
					y += mult
				else:		
					if k_ss >= 0:
						y += (valuation(k-k_ss,p)+e)*mult
					if k_ss < 0:
						y += mult
			NP += [(i,y)]

		if num==None:
			return NewtonPolygon(NP).slopes()
		else:
			return NewtonPolygon(NP).slopes()[0:num]

	def multiplicity(self,comp,i):
		"""Returns the total number of zeroes with multiplicity in the i-th coefficient
			on component comp"""
		return sum([self[comp][i][a][1] for a in range(len(self[comp][i]))])

	def wadic_slopes(self,comp,num=None):
		"""Returns the w-adic slopes of the mod p reduction of the ghost series on 
			component comp"""
		if num!=None:
			##HACKING HERE ABOUT NUMBER OF TERMS
			if self.num_coefs(comp)<num+10:
				self.compute_ghost_series(num_coefs=num+10)

		NP = [(a,self.multiplicity(comp,a)) for a in range(self.num_coefs(comp))]
		if num!=None:
			return NewtonPolygon(NP).slopes()[0:num]
		else:
			return NewtonPolygon(NP).slopes()

	def form_p2_modification(self):
		N = self.N
		new_ghost_0 = [x for x in self[0]]
		new_ghost = [new_ghost_0]
		max_i = len(new_ghost[0])
		dim_gaps = N*prod([1 + 1/ell for ell in ZZ(N).prime_factors()])
		chi = [c for c in DirichletGroup(8) if c.conductor() == 8 and c(-1) == 1][0]
		dim2 = dimension_cusp_forms(DirichletGroup(8*N)(c),2)
		max_k = floor((max_i - dim2)/dim_gaps + 2)
		mod_zeros = get_modified_zeros(N,max_k)
		for j in range(max_i):
			new_ghost[0][j] += mod_zeros[j]
		self.series=new_ghost

	#############################
	#SPECTRAL CURVE PICTURE STUFF
	#############################

	## i = index
	## v = valuation (non-integral)
	## comp = component
	## Return (r,s) where r = total ghost zeroes (with mult) > v
	## and s = sum v_p(w) where w is a ghost zero with v_p(w) < v. Thus
	## on the valuation v_p(w) = v we see v_p(a_i) = v*r + s.
	def zero_val_count(self,comp,i,v,central_wt=0):
		assert floor(v)!=v, "Use non-integral valuations"
		p = self.p
		ai = self[comp][i]
		if p == 2:
			e = 2  ## this is the extra valuation from changing from k to w
		else:
			e = 1
		r= sum([a[1] for a in ai if ZZ(a[0]-central_wt).valuation(p)+e > v])
		s= sum([(ZZ(a[0]-central_wt).valuation(p)+ e)*a[1] for a in ai if ZZ(a[0]-central_wt).valuation(p)+e < v])
		return (r,s)

	## v = valuation
	## Returns a list (r_i,s_i) such that (i,r_i*v+s_i) is the i-th Newton point of U_p in valuation v = v2(wt-central_wt)
	## added 5/15/16: ability to use a central_wt
	def newton_points_at_fixed_valuation(self,comp,v,max_i,central_wt=0):
		p = self.p
		return [self.zero_val_count(comp,i,v,central_wt) for i in range(max_i+1)]

	def global_halo_piece(self,comp,v,max_i,central_wt=0):
		p=self.p
		max_i = 2*max_i
		list = self.newton_points_at_fixed_valuation(comp,v,max_i,central_wt)
		left = NewtonPolygon([(i,list[i][0]*floor(v)+list[i][1]) for i in range(max_i+1)]).slopes()
		right = NewtonPolygon([(i,list[i][0]*ceil(v)+list[i][1]) for i in range(max_i+1)]).slopes()
		max_i = max_i/2
		r = floor(v)
		spectral = point((0,0))
		prev = [(r,0),(r+1,0)]
		linethickness = .1
		for i in range(max_i):
			linecolor = 'blue'
			linethickness = len([j for j in range(max_i) if left[j]==left[i] and right[j] == right[i]])
			spectral += halo_line((r,left[i]),(r+1,right[i]),linethickness,linecolor,prime=p)
			prev = [(floor(v),left[i]),(ceil(v),right[i])]
		return spectral

	## comp = component
	## sheets = number of sheets computed
	## max_v = maximal valuation of weights used
	## min_v (optional) = minimum valuation of weights used
	## central_wt (optional) = weight placed at center of weight space
	def global_halo(self,comp,sheets,max_v,min_v=0,central_wt = 0,xaxes=None,yaxes=None):
		p = self.p
		max_i = sheets
		assert floor(min_v)==min_v and floor(max_v)==max_v, "C'mon.  Use an integer for the endpoint"
		if 2*max_i + 10 > self.num_coefs(comp):
			g.compute_ghost_series(2*max_i+10)
		gh = self.global_halo_piece(comp,min_v+0.5,max_i,central_wt)
		for v in range(min_v+1,max_v):
			gh +=  self.global_halo_piece(comp,v+0.5,max_i,central_wt)
		if xaxes != None or yaxes != None:
			gh.axes_range(0,xaxes,0,yaxes)

		return gh






####Helper function for halos

def halo_line(pt1,pt2,linethickness,linecolor,prime):
	p = prime
	r = 10
	outline = 'dashed'
	colored = 'white'
	if p >= 3 or (p == 2 and pt2[0] >= 4):
		return line([pt1,pt2],thickness=linethickness,color=linecolor) + point(pt1,size=linethickness*r,color=colored,faceted = True,alpha=1,zorder=10) + point(pt2,size=linethickness*r,color=colored,faceted = True,alpha=1,zorder=10)
	if p == 2 and pt2[0] < 4:
		return line([pt1,pt2],thickness=linethickness,color=linecolor)


        

## Helper functions for p=2 modification

def mults(v):
	mult_free = list(set(v))
	mult_free.sort()
	ans = []
	for a in mult_free:
		ans.append([a,v.count(a)])
	return ans
 
## takes in a list [(a0,d0),(a1,d1),...] where
## 0 = a0  < a1 < ... and di > 0
## and returns the correct multiplicity pattern
## for the weight k = 2 points
## sage: wt2_slopes_to_mults([(0,7),(1/2,1),(4,3)])
## [0, 1, 2, 3, 3, 2, 1, 0, 0, 1, 1, 0]
 
def wt2_slopes_to_mults(mult_slope_list):
	mults = []
	ell = len(mult_slope_list)
	for i in range(ell):
		dimi = mult_slope_list[i][1]
		for j in range(1,dimi+1):
			if j < dimi/2:
				mults.append(j)
			else:
				mults.append(dimi+1 - (j+1))
	return mults
	

### returns a list [(-k,m)]_i where we want to make a modification
### to the ghost a_i by adding in a zero of mult m at the weight -5^k - 1.

def get_modified_zeros(N,max_k):
	# get even character of conductor 8
	chi = [c for c in DirichletGroup(8) if c.conductor() == 8 and c(-1) == 1][0]
	S = ModularSymbols(DirichletGroup(N*8)(chi),weight=2,sign=1).cuspidal_subspace()
	#S = CuspForms(DirichletGroup(N*8)(chi),2)
	# computes weight 2 slopes with cond. 8
	mult_slope_list = mults(S.hecke_polynomial(2).newton_slopes(2))
	# isolates just the fractional slopes in the classical space
	base_mults = wt2_slopes_to_mults([x for x in mult_slope_list if x[0] != 0 and x[0] != 1])
	if mult_slope_list[0][0] == 0:
		start_ind = mult_slope_list[0][1] + 1
	else:
		start_ind = 1
	gap = (S.dimension() + Gamma0(N).ncusps())
	mods = [[]]*(start_ind + len(base_mults) + gap*(max_k-1))
	for k in range(2,max_k+1):
		for j in range(len(base_mults)):
			if base_mults[j] != 0:
				mods[j+start_ind + gap*(k-2)] = [(-k,base_mults[j])]
				#mods[j+start_ind + gap*(k-2)] = [(-2,base_mults[j])]
	return mods

### fun fact: if N is odd then
### dim S_{k+1}(Gamma1(8N),e\pm) - dim S_k(Gamma1(8N),e\pm) = N\cdot prod_{\ell \dvd N} (1 + 1/\ell)

def form_p2_modification(N,ghost):
	new_ghost_0 = [x for x in ghost[0]]
	new_ghost = [new_ghost_0,[]]
	max_i = len(new_ghost[0])
	dim_gaps = N*prod([1 + 1/ell for ell in ZZ(N).prime_factors()])
	chi = [c for c in DirichletGroup(8) if c.conductor() == 8 and c(-1) == 1][0]
	dim2 = dimension_cusp_forms(DirichletGroup(8*N)(c),2)
	max_k = floor((max_i - dim2)/dim_gaps + 2)
	mod_zeros = get_modified_zeros(N,max_k)
	for j in range(max_i):
		new_ghost[0][j] += mod_zeros[j]
	return new_ghost
	
def spread_out(mult_list):
	new_list = []
	for x in mult_list:
		for j in range(x[1]):
			new_list.append(x[0])
	return new_list
	
def gather_fractional_slopes(slope_data):
	ret_data = []
	for y in slope_data:
		wt = y[0]
		twist_data = []
		for z in y[1]:
			twist = z[0]
			mult_list = z[1]
			new_mults = spread_out(mult_list)
			frac_data = [(x,new_mults.index(x)) for x in set(new_mults) if not x.is_integer()]
			twist_data.append([twist,frac_data])
		ret_data.append([wt,twist_data])
	return ret_data
