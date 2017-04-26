__kernel void mdh(  __global float *ax, __global float *ay, __global float *az,
					__global float *charge, __global float *size, 
					__global float *gx,	__global float *gy, __global float *gz,
					float pre1, float xkappa, __global float *val, int natoms, 
					__local float * shared  )
{
	   int igrid = get_global_id(0);
	   int iatom;
	   float v = 0.0f;
	   
	   float lgx = gx[igrid];
	   float lgy = gy[igrid];
	   float lgz = gz[igrid];
	   
	   for( iatom = 0; iatom < natoms; iatom++ )
	   {
			 float dx = lgx - ax[iatom];
			 float dy = lgy - ay[iatom];
			 float dz = lgz - az[iatom];

			 float dist = sqrt( dx * dx + dy * dy + dz * dz );
			 v += pre1 * ( charge[iatom] / dist )  * 
					exp( -xkappa * (dist - size[iatom])) / 
					(1.0f + xkappa * size[iatom]);
	   }

	   val[ igrid ] = v;
}


__kernel void plusone(	__global unsigned long* reactantSpeciesIndices,
						__global unsigned long int nReactants)
{
	   int global_id = get_global_id(0);

	   unsigned long index = reactantSpeciesIndices[global_id];
	   reactantSpeciesIndices[global_id]=index+1;

}



__kernel void computeTotalPropensity ( 
		__global unsigned long* species,
		__global unsigned long* reactantSpeciesIndices,
		__global unsigned long* reactantSpeciesStoichiometry,
		__global unsigned long* reactantSpeciesHOR,
		__global float* propensities)
{
	   int global_id = get_global_id(0);
		
	
}

__kernel void computeCumulativePropensities (
		__global float* propensities,
		__global float* cumulativePropensities,
		int nReactions)
{
	int global_id = get_global_id(0);
	
}


__kernel void  computeMuHatSigmaHat(
		__global float* propensities,
		__global unsigned long* reactantSpeciesStoichiometry,
		__global float* muHat,
		__global float* sigmaHat
		)
{
	   int global_id = get_global_id(0);
}

__kernel void  computeLeapSize(
		__global float* muHat,
		__global float* sigmaHat,
		__global float* xi,
		__global unsigned long* L1,
		__global unsigned long* L2,
		__global unsigned long* L)
{
	   int global_id = get_global_id(0);
}

__kernel void  updateSpecies(
		__global unsigned long* species,
		__global unsigned long* proposedSpecies,
		__global unsigned long* L,
		__global float* randomNumbers)
{
	   int global_id = get_global_id(0);
}

__kernel void  compute_xi(	__global unsigned long* reactantSpeciesHOR,
							__global unsigned long* reactantSpeciesHnu,
							__global unsigned long* reactantSpeciesStoichiometry,
							__global unsigned long* reactantSpecies,
							__global float* xi,
							int nSpecies,
							float epsilonInv	)
{
	int i = get_global_id(0);
	if (i<nSpecies)
	{
		float x=reactantSpecies[i];
		switch (reactantSpeciesHOR[i]) {
			case 2:
				switch (reactantSpeciesHnu[i]) {
					case 2:
						x/=(2+(1/(x-1)));
						break;
					case 1:
						x/=2.0f;
						break;
				}
				break;
			case 3:
				switch (reactantSpeciesHnu[i]) {
					case 1:
						x/=3;
						break;
					case 2:
						x/=((3+(3/(2*x-2))));
						break;
					case 3:
						x/=(3+(1/(x-1))+(2/(x-2)));
				}
		}
		
		if (x>epsilonInv)
		{
			xi[i]=1.0f;
		}
		else
		{
			xi[i]=x/epsilonInv;
		}

	}
	
}

__kernel void computeLj(
		__global unsigned long* species,
		__global unsigned long* reactantSpeciesIndices,
		__global unsigned long* Lj)
{
	   int global_id = get_global_id(0);
}