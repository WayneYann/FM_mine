Double TSoot::Nucleation( int i, Double temp, Double *pahMoments )
{
	switch( i ) {
		case 0: 
			return 0.5 * GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[0];
		case 1:
			return GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[1];
		case 2:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[2] 
								+ pahMoments[1] * pahMoments[1] );
		case 3:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[3] 
								+ 3.0 * pahMoments[1] * pahMoments[2] );
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
}

Double TSoot::NucleationNow( int i, Double temp, Double *Y, Double rho, Double *molarMass )
{	
	switch( i ) {
		case 0: 
			return 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		case 1:
			return 18.0 * 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		case 2:
			return 324.0 * 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		case 3:
			return 5832.0 * 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
}

//Mueller (9/7/07)
#ifdef VS
Double TSoot::NucleationNew(int k, Double temp, Double *pahMoments)
{
	switch (k)
	{
		case 0:	//M(0,0)
			return 0.5 * GetPhi(temp, pahMoments, pahMoments, 1.0/2.0); //Coalescence
		case 1:	//M(1,0)
			return GetPhi(1, 0, 0, 0, temp, 3.0, pahMoments, pahMoments); //Coalescence
		case 2:	//M(0,1)
			return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/6.0); //Coalescence
		case 3:	//M(-1/2,0)
			return 0.5 * GetPhi(temp, pahMoments, pahMoments, 0); //Coalescence
		case 4:	//M(0,2)
			return 0.5 * GetPhi(temp, pahMoments, pahMoments, 11.0/6.0); //Coalescence
		case 5:	//M(7/6,1)
			return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/3.0); //Coalescence
		default:
			cerr << "#error: nucleation not yet implemented for " << k << " moments" << endl;
			exit(2);
	}

	return 0.0;
}
#endif

//Mueller (7/19/07)
Double TSoot::NucleationNew(int i, int j, Double temp, Double *pahMoments)
{
	return 0.5 * GetPhi(temp, pahMoments, pahMoments, (double(i)/6.0)+(2.0/3.0)*(double(j)/6.0)+(1.0/2.0));

	if (i==0 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 1.0/2.0);
	else if (i==6 && j==0)
		return GetPhi(1, 0, 0, 0, temp, 3.0, pahMoments, pahMoments);
	else if (i==0 && j==6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/6.0);
	else if (i==0 && j==3)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 5.0/6.0);
	else if (i==2 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 5.0/6.0);
	else if (i==4 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/6.0);
	else if (i==2 && j==3)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/6.0);
	else if (i==6 && j==6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 13.0/6.0);
	else if (i==12 && j==0)
		return GetPhi(2, 0, 0, 0, temp, 3.0, pahMoments, pahMoments) + GetPhi(1, 0, 1, 0, temp, 3.0, pahMoments, pahMoments);
/*	else if (i==-12 && j==18)
		return 0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, pahMoments, pahMoments);*/
	else if (i==11 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/3.0);
	else if (i==7 && j==6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 7.0/3.0);
	else if (i==0 && j==18)
		return GetPhi(2, 0, 0, 0, temp, 3.0, pahMoments, pahMoments) + GetPhi(1, 0, 1, 0, temp, 3.0, pahMoments, pahMoments);
	else if (i==-3 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 0);
	else if (i==13 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 8.0/3.0);
	else if (i==7 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 5.0/3.0);
	else if (i==4 && j==3)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 3.0/2.0);
//		return GetPhi(1, 0, 0, 0, temp, 3.0, pahMoments, pahMoments);
	else if (i==0 && j==-6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -1.0/6.0);
	else if (i==6 && j==-6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 5.0/6.0);
	else if (i==12 && j==-6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 11.0/6.0);
	else if (i==-6 && j==-6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -7.0/6.0);
	else if (i==-6 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -1.0/2.0);
	else if (i==-6 && j==6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 1.0/6.0);
	else if (i==-6 && j==12)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 5.0/6.0);
	else if (i==-12 && j==-6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -13.0/6.0);
	else if (i==-12 && j==0)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -3.0/2.0);
	else if (i==-12 && j==6)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -5.0/6.0);
	else if (i==-12 && j==12)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, -1.0/6.0);
	else if (i==-12 && j==18)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 1.0/2.0);
	else if (i==0 && j==12)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 11.0/6.0);
/*	else if (i==2 && j==0)
		return GetPhi(2, 0, 0, 0, temp, 3.0, pahMoments, pahMoments) + GetPhi(1, 0, 1, 0, temp, 3.0, pahMoments, pahMoments);
	else if (i==0 && j==2)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 11.0/6.0);
	else if (i==1 && j==1)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 13.0/6.0);
	else if (i==3 && j==0)
		return GetPhi(3, 0, 0, 0, temp, 3.0, pahMoments, pahMoments) + 3.0 * GetPhi(2, 0, 1, 0, temp, 3.0, pahMoments, pahMoments);
	else if (i==0 && j==3)
		return GetPhi(2, 0, 0, 0, temp, 3.0, pahMoments, pahMoments) + GetPhi(1, 0, 1, 0, temp, 3.0, pahMoments, pahMoments);
	else if (i==2 && j==1)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 19.0/6.0);
	else if (i==1 && j==2)
		return 0.5 * GetPhi(temp, pahMoments, pahMoments, 17.0/6.0);*/
	else
	{
		cerr << "#error: no need to compute source term for nucleation for M_" << i << "," << j << endl;
		exit(2);
	}

	return 0.0;
}

#ifndef VS
Double TSoot::NucleationNew( int i, Double temp, Double *pahMoments )
{	
#ifdef NEWPOLYNUCCOND
	switch( i ) {
		case 0: 
			return 0.5 * GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[0];
		case 1:
			return GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[1];
		case 2:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[2] 
								+ pahMoments[1] * pahMoments[1] );
		case 3:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[3] 
								+ 3.0 * pahMoments[1] * pahMoments[2] );
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
#else
	switch( i ) {
#	ifdef NOPAH
#		ifdef PAHFROMA4
		case 0: 
			return 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
		case 1:
			return 18.0 * 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
		case 2:
			return 324.0 * 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
		case 3:
			return 5832.0 * 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
#		else
			fprintf( stderr, "#error in TSoot::NucleationNew: not yet implemented\n" );
			exit( 2 );
#		endif
#	else
		case 0: 
		  return 0.5 * GetPhiPAH( 0, 0, temp, pahMoments );//changed by GB
		case 1:
		  return GetPhiPAH( 0, 1, temp, pahMoments );//changed by GB 
		case 2:
			return GetPhi( 0, 2, temp, pahMoments ) + GetPhi( 1, 1, temp, pahMoments );
		case 3:
			return GetPhi( 0, 3, temp, pahMoments ) + 3.0 * GetPhi( 1, 2, temp, pahMoments );
#	endif
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
#endif
}
#endif

//Mueller (8/9/07)
Double TSoot::FractionalMoment2(Double Index1, Double Index2, Double *moments)
{
	//M(0,0), M(1,0), M(0,1)
	Double a = pow(MAX(moments[1], 1.0e-30), 1.0) * pow(MAX(moments[0], 1.0e-30), -1.0);
	Double b = pow(MAX(moments[2], 1.0e-30), 1.0) * pow(MAX(moments[0], 1.0e-30), -1.0);
	Double c = MAX(moments[0], 1.0e-30);

	Double aexp = pow(a, Index1);
	Double bexp = pow(b, Index2);
	Double cexp = c;

	return aexp * bexp * cexp;
}

//Mueller (9/18/07)
Double TSoot::FractionalMoment3(Double Index1, Double Index2, Double *moments)
{
	//Assumed distribution is two delta distributions.
	//Functional form of moments is: Mxy=a*exp(bx+cy)+d*exp(ex+fy).
	//Coefficients found numerically.
	//First three moments must be M(0,0), M(1,0), M(0,1).

	double * c = new double[fNSootMoments];

	double variance = (MAX(moments[3],1.0e-30)/MAX(moments[0],1.0e-30)) - pow(MAX(moments[1],1.0e-30)/MAX(moments[0],1.0e-30),double(Geti(3))/6.0) * pow(MAX(moments[2],1.0e-30)/MAX(moments[0],1.0e-30),double(Getj(3))/6.0);
	double normvar = variance / (MAX(moments[3],1.0e-30)/MAX(moments[0],1.0e-30));

	if (fabs(normvar) < 1.0e-06)
	{
		c[0] = 0.5 * MAX(moments[0],1.0e-30);
		c[1] = log(MAX(moments[1],1.0e-30) / MAX(moments[0],1.0e-30));
		c[2] = log(MAX(moments[2],1.0e-30) / MAX(moments[0],1.0e-30));
		c[3] = 0.5 * MAX(moments[0],1.0e-30);
		c[4] = log(MAX(moments[1],1.0e-30) / MAX(moments[0],1.0e-30));
		c[5] = log(MAX(moments[2],1.0e-30) / MAX(moments[0],1.0e-30));
	}
	else
	{
		double * c_curr = new double[fNSootMoments];
		double * c_next = new double[fNSootMoments];

		c_curr[0] = (1.0/2.0) * moments[0];
		c_curr[1] = log(18.0);
		c_curr[2] = log(pow(18.0, 2.0/3.0));
		c_curr[3] = (1.0/2.0) * moments[0];
		c_curr[4] = log(2.0 * (moments[1] / moments[0]) - 18.0);
		c_curr[5] = log(2.0 * (moments[2] / moments[0]) - pow(18.0, 2.0/3.0));

		double * f = new double[fNSootMoments];
		double * mdeltax = new double[fNSootMoments];
	
		double ** J = new double*[fNSootMoments];
		for (int i = 0; i < fNSootMoments; i++)
			J[i] = new double[fNSootMoments];

		bool converge;

		do
		{
			FuncEval(f, c_curr, moments);

			FillJac(J, c_curr, moments);
	
			Solve(J, mdeltax, f);

			Subtract(c_next, c_curr, mdeltax);

			converge = Converge(c_next, c_curr);

			Assign(c_curr, c_next);
		}
		while (!converge);

		Assign(c, c_curr);
	
		delete [] f;

		delete [] c_curr;
		delete [] c_next;

		for (int i = 0; i < fNSootMoments; i++)
		{
			delete [] J[i];
		}
		delete [] J;

		delete [] mdeltax;
	}

	double fracmom = c[0]*exp(c[1]*Index1+c[2]*Index2)+c[3]*exp(c[4]*Index1+c[5]*Index2);

/*	cerr << "a: " << c[0] << endl;
	cerr << "b: " << c[1] << endl;
	cerr << "c: " << c[2] << endl;
	cerr << "d: " << c[3] << endl;
	cerr << "e: " << c[4] << endl;
	cerr << "f: " << c[5] << endl;
	cin.get();*/

/*	if (isnan(fracmom))
	{
		c[0] = c[3] = 0.5 * moments[0];
		c[1] = c[4] = log(MAX(moments[1], 1.0e-30) / MAX(moments[0], 1.0e-30));
		c[2] = c[5] = log(MAX(moments[2], 1.0e-30) / MAX(moments[0], 1.0e-30));

		fracmom = c[0]*exp(c[1]*Index1+c[2]*Index2)+c[3]*exp(c[4]*Index1+c[5]*Index);
	}*/

	delete [] c;
	
	return fracmom;
}

void TSoot::FuncEval(double * f, double * c, Double * moments)
{
	for (int i = 0; i < fNSootMoments; i++)
	{
		f[i] = c[0]*exp(c[1]*double(Geti(i))/6.0+c[2]*double(Getj(i))/6.0)+c[3]*exp(c[4]*double(Geti(i))/6.0+c[5]*double(Getj(i))/6.0)-moments[i];
	}
}

void TSoot::FillJac(double ** J, double * c, Double * moments)
{
	for (int i = 0; i < fNSootMoments; i++)
	{
		for (int j = 0; j < fNSootMoments; j++)
		{
			J[i][j] = DerivEval(i, j, c, moments);
		}
	}
}

double TSoot::DerivEval(int eq, int var, double * c, Double * moments)
{
	switch (var)
	{
		case 0:	return exp(c[1]*double(Geti(eq))/6.0+c[2]*double(Getj(eq))/6.0);
		case 1:	return c[0] * double(Geti(eq))/6.0 * exp(c[1]*double(Geti(eq))/6.0+c[2]*double(Getj(eq))/6.0);
		case 2:	return c[0] * double(Getj(eq))/6.0 * exp(c[1]*double(Geti(eq))/6.0+c[2]*double(Getj(eq))/6.0);
		case 3: return exp(c[4]*double(Geti(eq))/6.0+c[5]*double(Getj(eq))/6.0);
		case 4:	return c[3] * double(Geti(eq))/6.0 * exp(c[4]*double(Geti(eq))/6.0+c[5]*double(Getj(eq))/6.0);
		case 5:	return c[3] * double(Getj(eq))/6.0 * exp(c[4]*double(Geti(eq))/6.0+c[5]*double(Getj(eq))/6.0);
	}

	return 0.0;
}

void TSoot::Solve(double ** A, double * x, double * b)
{
	double ** Aug = new double*[fNSootMoments];
	for (int i = 0; i < fNSootMoments; i++)
		Aug[i] = new double[fNSootMoments+1];

	for (int i = 0; i < fNSootMoments; i++)
	{
		for (int j = 0; j < fNSootMoments; j++)
		{
			Aug[i][j] = A[i][j];
		}
		Aug[i][fNSootMoments] = b[i];
	}

	int i = 0;

	do
	{
		double temp = Aug[i][i];

		for (int j = i; j < fNSootMoments+1; j++)
		{
			Aug[i][j] /= temp;
		}

		for (int j = 0; j < fNSootMoments; j++)
		{
			if (j != i)
			{
				temp = Aug[j][i];
				for (int k = i; k < fNSootMoments+1; k++)
					Aug[j][k] -= temp*Aug[i][k];
			}
		}

		i++;
	}
	while (i < fNSootMoments);

	for (int i = 0; i < fNSootMoments; i++)
		x[i] = Aug[i][fNSootMoments];

	for (int i = 0; i < fNSootMoments; i++)
		delete [] Aug[i];
	delete [] Aug;
}

void TSoot::Subtract(double * diff, double * min, double * sub)
{
	for (int i = 0; i < fNSootMoments; i++)
		diff[i] = min[i] - sub[i];
}

bool TSoot::Converge(double * next, double * curr)
{
	bool converge = true;
/*	bool smalldelta = false;
	bool largedelta = false;

	if (fabs(next[2]) < 1.0e-8 && fabs(curr[2]) < 1.0e-8)
	{
		smalldelta = true;
		next[2] = 0.0;
	}
	
	if (fabs(next[5]) < 1.0e-8 && fabs(curr[5]) < 1.0e-8)
	{
		largedelta = true;
		next[5] = 0.0;
	}*/

	for (int i = 0; i < fNSootMoments; i++)
	{
//		if (smalldelta && i==2)
//			continue;
//		if (largedelta && i==5)
//			continue;
		if (fabs(next[i] - curr[i]) / fabs(curr[i]) > 1.0e-8)
			converge = false;
	}

	return converge;	
}

void TSoot::Assign(double * left, double * right)
{
	for (int i = 0; i < fNSootMoments; i++)
		left[i] = right[i];
}

#ifdef MUELLER
//Mueller (8/14/07)
Double TSoot::FractionalMoment3(Double Index1, Double Index2, Double *moments)
{
	//M(0,0), M(0,1), M(1,0), M(0,2), M(1,1), M(2,0)
/*	Double a = pow(MAX(moments[5], 1.0e-30), 1.0/2.0) * pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[2], 1.0e-30), -1.0);
	Double b = pow(MAX(moments[3], 1.0e-30), 1.0/2.0) * pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[1], 1.0e-30), -1.0);
	Double c = pow(MAX(moments[2], 1.0e-30), 2.0) * pow(MAX(moments[5], 1.0e-30), -1.0/2.0) * pow(MAX(moments[0], 1.0e-30), -3.0/2.0);
	Double d = pow(MAX(moments[1], 1.0e-30), 2.0) * pow(MAX(moments[3], 1.0e-30), -1.0/2.0) * pow(MAX(moments[0], 1.0e-30), -3.0/2.0);
	Double e = MAX(moments[0], 1.0e-30);
	Double f = pow(MAX(moments[4], 1.0e-30), 1.0) * pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[2], 1.0e-30), -1.0) * pow(MAX(moments[1], 1.0e-30), -1.0);*/

/*	Double a = pow(moments[5], 1.0/2.0) * pow(moments[0], 1.0/2.0) * pow(moments[2], -1.0);
	Double b = pow(moments[3], 1.0/2.0) * pow(moments[0], 1.0/2.0) * pow(moments[1], -1.0);
	Double c = pow(moments[2], 2.0) * pow(moments[5], -1.0/2.0) * pow(moments[0], -3.0/2.0);
	Double d = pow(moments[1], 2.0) * pow(moments[3], -1.0/2.0) * pow(moments[0], -3.0/2.0);
	Double e = moments[0];
	Double f = pow(moments[4], 1.0) * pow(moments[0], 1.0) * pow(moments[2], -1.0) * pow(moments[1], -1.0);*/

/*if (moments[0] < 1.0e-20 || moments[1] < 1.0e-20 || moments[2] < 1.0e-20 || moments[3] < 1.0e-20 || moments[4] < 1.0e-20 || moments[5] < 1.0e-20)
{
//	cerr << "FracMom clipping...\n";
	return 1.0e-20;
}*/

/*	Double a = pow(moments[2], 3.0) * pow(moments[4], -4.5) * pow(moments[0], 1.5);
	Double b = pow(moments[4], 4.5) * pow(moments[2], -2.5) * pow(moments[0], -2.5);
	Double c = pow(moments[5], 6.0) * pow(moments[3], -6.0) * pow(moments[2], 2.0) * pow(moments[4], -6.0) * pow(moments[0], 4.0);
	Double d = pow(moments[3], 4.0) * pow(moments[1], -1.0) * pow(moments[0], -3.0);
	Double e = pow(moments[1], 2.0) * pow(moments[3], -4.0) * pow(moments[0], 2.0);
	Double f = moments[0];*/

/*	Double aexp = pow(a, Index1*Index1);
	Double bexp = pow(b, Index1);
	Double cexp = pow(c, Index1*Index2);
	Double dexp = pow(d, Index2);
	Double eexp = pow(e, Index2*Index2);
	Double fexp = f;

	return aexp * bexp * cexp * dexp * eexp * fexp;*/

/*	Double * temp;
	temp = new Double[fNSootMoments];
//cerr << "temp declared and allocated...\n";
	for (int i = 0; i < fNSootMoments; ++i)
	{
		temp[i] = log(MAX(moments[i], 1.0e-30));
//		cerr << "Moment (" << i << "): " << moments[i] << endl;
//		cerr << "Log of Moment (" << i << "): " << temp[i] << endl;
	}
//cin.get();
	Double a = 3.0 * temp[2] - 4.5 * temp[4] + 1.5 * temp[0];
	Double b = 4.5 * temp[4] - 2.0 * temp[2] - 2.5 * temp[0];
	Double c = 6.0 * temp[5] - 6.0 * temp[3] + 2.0 * temp[2] - 6.0 * temp[4] + 4.0 * temp[0];
	Double d = 4.0 * temp[3] - 1.0 * temp[1] - 3.0 * temp[0];
	Double e = 2.0 * temp[1] - 4.0 * temp[3] + 2.0 * temp[0];
	Double f = temp[0];

	Double aexp = a * Index1 * Index1;
	Double bexp = b * Index1;
	Double cexp = c * Index1 * Index2;
	Double dexp = d * Index2;
	Double eexp = e * Index2 * Index2;
	Double fexp = f;

	delete temp;

	return exp(aexp + bexp + cexp + dexp + eexp + fexp);*/
	
/*	Double a = pow(MAX(moments[2], 1.0e-20), 3.0) * pow(MAX(moments[4], 1.0e-20), -4.5) * pow(MAX(moments[0], 1.0e-20), 1.5);
	Double b = pow(MAX(moments[4], 1.0e-20), 4.5) * pow(MAX(moments[2], 1.0e-20), -2.0) * pow(MAX(moments[0], 1.0e-20), -2.5);
	Double c = pow(MAX(moments[5], 1.0e-20), 6.0) * pow(MAX(moments[3], 1.0e-20), -6.0) * pow(MAX(moments[2], 1.0e-20), 2.0) * pow(MAX(moments[4], 1.0e-20), -6.0) * pow(MAX(moments[0], 1.0e-20), 4.0);
	Double d = pow(MAX(moments[3], 1.0e-20), 4.0) * pow(MAX(moments[1], 1.0e-20), -1.0) * pow(MAX(moments[0], 1.0e-20), -3.0);
	Double e = pow(MAX(moments[1], 1.0e-20), 2.0) * pow(MAX(moments[3], 1.0e-20), -4.0) * pow(MAX(moments[0], 1.0e-20), 2.0);
	Double f = MAX(moments[0], 1.0e-20);*/

	//M(0,0), M(1,0), M(0,1), M(0,1/2), M(2/3,0), M(1,1)
/*	Double a = pow(MAX(moments[2], 1.0e-30), 3.0) * pow(MAX(moments[4], 1.0e-30), -4.5) * pow(MAX(moments[0], 1.0e-30), 1.5);
	Double b = pow(MAX(moments[4], 1.0e-30), 4.5) * pow(MAX(moments[2], 1.0e-30), -2.0) * pow(MAX(moments[0], 1.0e-30), -2.5);
	Double c = pow(MAX(moments[5], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), -1.0) * pow(MAX(moments[2], 1.0e-30), -1.0) * pow(MAX(moments[0], 1.0e-30), 1.0);
	Double d = pow(MAX(moments[3], 1.0e-30), 4.0) * pow(MAX(moments[1], 1.0e-30), -1.0) * pow(MAX(moments[0], 1.0e-30), -3.0);
	Double e = pow(MAX(moments[1], 1.0e-30), 2.0) * pow(MAX(moments[3], 1.0e-30), -4.0) * pow(MAX(moments[0], 1.0e-30), 2.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(2,0), M(0,1), M(1,1), M(-2,3)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[1], 1.0e-30), -1.0) * pow(MAX(moments[2], 1.0e-30), 1.0/2.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[1], 1.0e-30), 2.0) * pow(MAX(moments[2], 1.0e-30), -1.0/2.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), -1.0) * pow(MAX(moments[3], 1.0e-30), -1.0) * pow(MAX(moments[4], 1.0e-30), 1.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[1], 1.0e-30), -1.0/3.0) * pow(MAX(moments[2], 1.0e-30), 1.0/2.0) * pow(MAX(moments[3], 1.0e-30), 5.0/2.0) * pow(MAX(moments[4], 1.0e-30), -1.0) * pow(MAX(moments[5], 1.0e-30), -1.0/6.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[1], 1.0e-30), 1.0/3.0) * pow(MAX(moments[2], 1.0e-30), -1.0/2.0) * pow(MAX(moments[3], 1.0e-30), -3.0/2.0) * pow(MAX(moments[4], 1.0e-30), 1.0) * pow(MAX(moments[5], 1.0e-30), 1.0/6.0);
	Double f = pow(MAX(moments[0], 1.0e-30), 1.0);*/
	
	//M(0,0), M(1,0), M(11/6,0), M(0,1), M(7/6,1), M(-2,3)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 6.0/11.0) * pow(MAX(moments[1], 1.0e-30), -6.0/5.0) * pow(MAX(moments[2], 1.0e-30), 36.0/55.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -17.0/11.0) * pow(MAX(moments[1], 1.0e-30), 11.0/5.0) * pow(MAX(moments[2], 1.0e-30), -36.0/55.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 10.0/11.0) * pow(MAX(moments[1], 1.0e-30), -4.0/5.0) * pow(MAX(moments[2], 1.0e-30), -6.0/55.0) * pow(MAX(moments[3], 1.0e-30), -6.0/7.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -15.0/11.0) * pow(MAX(moments[1], 1.0e-30), -11.0/15.0) * pow(MAX(moments[2], 1.0e-30), 42.0/55.0) * pow(MAX(moments[3], 1.0e-30), 33.0/14.0) * pow(MAX(moments[4], 1.0e-30), -6.0/7.0) * pow(MAX(moments[5], 1.0e-30), -1.0/6.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 4.0/11.0) * pow(MAX(moments[1], 1.0e-30), 11.0/15.0) * pow(MAX(moments[2], 1.0e-30), -42.0/55.0) * pow(MAX(moments[3], 1.0e-30), -19.0/14.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0) * pow(MAX(moments[5], 1.0e-30), 1.0/6.0);
	Double f = pow(MAX(moments[0], 1.0e-30), 1.0);*/

	//M(0,0), M(1,0), M(0,1), M(2,0), M(1,1), M(0,3)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[1], 1.0e-30), -1.0) * pow(MAX(moments[3], 1.0e-30), 1.0/2.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[1], 1.0e-30), 2.0) * pow(MAX(moments[3], 1.0e-30), -1.0/2.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), -1.0) * pow(MAX(moments[2], 1.0e-30), -1.0) * pow(MAX(moments[4], 1.0e-30), 1.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -4.0/3.0) * pow(MAX(moments[2], 1.0e-30), 3.0/2.0) * pow(MAX(moments[5], 1.0e-30), -1.0/6.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 1.0/3.0) * pow(MAX(moments[2], 1.0e-30), -1.0/2.0) * pow(MAX(moments[5], 1.0e-30), 1.0/6.0);
	Double f = pow(MAX(moments[0], 1.0e-30), 1.0);*/

	//M(0,0), M(1,0), M(0,1), M(-1/2,0), M(7/6,1), M(1/3,1/2)
/*	Double a = pow(MAX(moments[0], 1.0e-30), -2.0) * pow(MAX(moments[1], 1.0e-30), 2.0/3.0) * pow(MAX(moments[3], 1.0e-30), 4.0/3.0);
	Double b = pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), 1.0/3.0) * pow(MAX(moments[3], 1.0e-30), -4.0/3.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 4.0/3.0) * pow(MAX(moments[1], 1.0e-30), -10.0/9.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -2.0/9.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -13.0/3.0) * pow(MAX(moments[2], 1.0e-30), -3.0/7.0) * pow(MAX(moments[3], 1.0e-30), 4.0/3.0) * pow(MAX(moments[4], 1.0e-30), -4.0/7.0) * pow(MAX(moments[5], 1.0e-30), 4.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 10.0/3.0) * pow(MAX(moments[2], 1.0e-30), 10.0/7.0) * pow(MAX(moments[3], 1.0e-30), -4.0/3.0) * pow(MAX(moments[4], 1.0e-30), 4.0/7.0) * pow(MAX(moments[5], 1.0e-30), -4.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(13/6,0), M(7/6,1), M(1/3,1/2)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 6.0/13.0) * pow(MAX(moments[1], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), 36.0/91.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -19.0/13.0) * pow(MAX(moments[1], 1.0e-30), 13.0/7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/91.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 12.0/13.0) * pow(MAX(moments[1], 1.0e-30), -6.0/7.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -6.0/91.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -73.0/39.0) * pow(MAX(moments[1], 1.0e-30), -32.0/21.0) * pow(MAX(moments[2], 1.0e-30), -3.0/7.0) * pow(MAX(moments[3], 1.0e-30), 36.0/91.0) * pow(MAX(moments[4], 1.0e-30), -4.0/7.0) * pow(MAX(moments[5], 1.0e-30), 4.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 34.0/39.0) * pow(MAX(moments[1], 1.0e-30), 32.0/21.0) * pow(MAX(moments[2], 1.0e-30), 10.0/7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/91.0) * pow(MAX(moments[4], 1.0e-30), 4.0/7.0) * pow(MAX(moments[5], 1.0e-30), -4.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(7/6,0), M(7/6,1), M(1/3,1/2)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 6.0/7.0) * pow(MAX(moments[1], 1.0e-30), -6.0) * pow(MAX(moments[3], 1.0e-30), 36.0/7.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -13.0/7.0) * pow(MAX(moments[1], 1.0e-30), 7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/7.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 6.0/7.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -6.0/7.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -31.0/21.0) * pow(MAX(moments[1], 1.0e-30), -20.0/3.0) * pow(MAX(moments[2], 1.0e-30), -3.0/7.0) * pow(MAX(moments[3], 1.0e-30), 36.0/7.0) * pow(MAX(moments[4], 1.0e-30), -4.0/7.0) * pow(MAX(moments[5], 1.0e-30), 4.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 10.0/21.0) * pow(MAX(moments[1], 1.0e-30), 20.0/3.0) * pow(MAX(moments[2], 1.0e-30), 10.0/7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/7.0) * pow(MAX(moments[4], 1.0e-30), 4.0/7.0) * pow(MAX(moments[5], 1.0e-30), -4.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(-1/2,0), M(7/6,1), M(2/3,1/2)
/*	Double a = pow(MAX(moments[0], 1.0e-30), -2.0) * pow(MAX(moments[1], 1.0e-30), 2.0/3.0) * pow(MAX(moments[3], 1.0e-30), 4.0/3.0);
	Double b = pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), 1.0/3.0) * pow(MAX(moments[3], 1.0e-30), -4.0/3.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 4.0/3.0) * pow(MAX(moments[1], 1.0e-30), -10.0/9.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -2.0/9.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -35.0/9.0) *pow(MAX(moments[1], 1.0e-30), -16.0/27.0) * pow(MAX(moments[2], 1.0e-30), 1.0/7.0) * pow(MAX(moments[3], 1.0e-30), 40.0/27.0) * pow(MAX(moments[4], 1.0e-30), -8.0/7.0) * pow(MAX(moments[5], 1.0e-30), 4.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 26.0/9.0) *pow(MAX(moments[1], 1.0e-30), 16.0/27.0) * pow(MAX(moments[2], 1.0e-30), 6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -40.0/27.0) * pow(MAX(moments[4], 1.0e-30), 8.0/7.0) * pow(MAX(moments[5], 1.0e-30), -4.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(13/6,0), M(7/6,1), M(2/3,1/2)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 6.0/13.0) * pow(MAX(moments[1], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), 36.0/91.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -19.0/13.0) * pow(MAX(moments[1], 1.0e-30), 13.0/7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/91.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 12.0/13.0) * pow(MAX(moments[1], 1.0e-30), -6.0/7.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -6.0/91.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -15.0/13.0) * pow(MAX(moments[1], 1.0e-30), -16.0/7.0) * pow(MAX(moments[2], 1.0e-30), 1.0/7.0) * pow(MAX(moments[3], 1.0e-30), 40.0/91.0) * pow(MAX(moments[4], 1.0e-30), -8.0/7.0) * pow(MAX(moments[5], 1.0e-30), 4.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 2.0/13.0) * pow(MAX(moments[1], 1.0e-30), 16.0/7.0) * pow(MAX(moments[2], 1.0e-30), 6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -40.0/91.0) * pow(MAX(moments[4], 1.0e-30), 8.0/7.0) * pow(MAX(moments[5], 1.0e-30), -4.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(7/6,0), M(7/6,1), M(1/3,1/2)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 6.0/7.0) * pow(MAX(moments[1], 1.0e-30), -6.0) * pow(MAX(moments[3], 1.0e-30), 36.0/7.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -13.0/7.0) * pow(MAX(moments[1], 1.0e-30), 7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/7.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 6.0/7.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -6.0/7.0) * pow(MAX(moments[4], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -5.0/7.0) * pow(MAX(moments[1], 1.0e-30), -8.0) * pow(MAX(moments[2], 1.0e-30), 1.0/7.0) * pow(MAX(moments[3], 1.0e-30), 40.0/7.0) * pow(MAX(moments[4], 1.0e-30), -8.0/7.0) * pow(MAX(moments[5], 1.0e-30), 4.0);
	Double e = pow(MAX(moments[0], 1.0e-30), -2.0/7.0) * pow(MAX(moments[1], 1.0e-30), 8.0) * pow(MAX(moments[2], 1.0e-30), 6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -40.0/7.0) * pow(MAX(moments[4], 1.0e-30), 8.0/7.0) * pow(MAX(moments[5], 1.0e-30), -4.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(13/6,0), M(0,2), M(7/6,1)
/*	Double a = pow(MAX(moments[0], 1.0e-30), 6.0/13.0) * pow(MAX(moments[1], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), 36.0/91.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -19.0/13.0) * pow(MAX(moments[1], 1.0e-30), 13.0/7.0) * pow(MAX(moments[3], 1.0e-30), -36.0/91.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 12.0/13.0) * pow(MAX(moments[1], 1.0e-30), -6.0/7.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -6.0/91.0) * pow(MAX(moments[5], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[2], 1.0e-30), 2.0) * pow(MAX(moments[4], 1.0e-30), -1.0/2.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[2], 1.0e-30), -1.0) * pow(MAX(moments[4], 1.0e-30), 1.0/2.0);
	Double f = MAX(moments[0], 1.0e-30);*/

	//M(0,0), M(1,0), M(0,1), M(-1/2,0), M(0,2), M(7/6,1)
/*	Double a = pow(MAX(moments[0], 1.0e-30), -2.0) * pow(MAX(moments[1], 1.0e-30), 2.0/3.0) * pow(MAX(moments[3], 1.0e-30), 4.0/3.0);
	Double b = pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), 1.0/3.0) * pow(MAX(moments[3], 1.0e-30), -4.0/3.0);
	Double c = pow(MAX(moments[0], 1.0e-30), 4.0/3.0) * pow(MAX(moments[1], 1.0e-30), -10.0/9.0) * pow(MAX(moments[2], 1.0e-30), -6.0/7.0) * pow(MAX(moments[3], 1.0e-30), -2.0/9.0) * pow(MAX(moments[5], 1.0e-30), 6.0/7.0);
	Double d = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[2], 1.0e-30), 2.0) * pow(MAX(moments[4], 1.0e-30), -1.0/2.0);
	Double e = pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[2], 1.0e-30), -1.0) * pow(MAX(moments[4], 1.0e-30), 1.0/2.0);
	Double f = MAX(moments[0], 1.0e-30);*/
	
/*	Double * temp;
	temp = new Double[fNSootMoments];
//cerr << "temp declared and allocated...\n";
	for (int i = 0; i < fNSootMoments; ++i)
	{
		if (moments[0] == 0 && moments[i] != 0)
		{
			cerr << "Houston, we have a problem...\n";
			cin.get();
		}
		if (moments[0] == 0)
			temp[i] = 1;
		else
		{
			temp[i] = moments[i] / moments[0];
//			cerr << "Temp(" << i << "): " << temp[i] << endl;
		}
//		cerr << "Moment (" << i << "): " << moments[i] << endl;
//		cerr << "Log of Moment (" << i << "): " << temp[i] << endl;
	}/*

/*	Double a = pow(temp[2], 3.0) * pow(temp[4], -4.5) * pow(temp[0], 1.5);
	Double b = pow(temp[4], 4.5) * pow(temp[2], -2.0) * pow(temp[0], -2.5);
	Double c = pow(temp[5], 6.0) * pow(temp[3], -6.0) * pow(temp[2], 2.0) * pow(temp[4], -6.0) * pow(temp[0], 4.0);
	Double d = pow(temp[3], 4.0) * pow(temp[1], -1.0) * pow(temp[0], -3.0);
	Double e = pow(temp[1], 2.0) * pow(temp[3], -4.0) * pow(temp[0], 2.0);
	Double f = temp[0];*/

/*	Double a = pow(temp[2], 3.0) * pow(temp[4], -4.5) * pow(temp[0], 1.5);
	Double b = pow(temp[4], 4.5) * pow(temp[2], -2.0) * pow(temp[0], -2.5);
	Double c = pow(temp[0], 1.0) * pow(temp[2], -1.0) * pow(temp[1], -1.0);
	Double d = pow(temp[3], 4.0) * pow(temp[1], -1.0) * pow(temp[0], -3.0);
	Double e = pow(temp[1], 2.0) * pow(temp[3], -4.0) * pow(temp[0], 2.0);
	Double f = temp[0];*/

/*	Double aexp = pow(a, Index1*Index1);
	Double bexp = pow(b, Index1);
	Double cexp = pow(c, Index1*Index2);
	Double dexp = pow(d, Index2);
	Double eexp = pow(e, Index2*Index2);
	Double fexp = f;*/

	//Two delta functions with M(0,0), M(1,0), M(0,1), M(2,0), M(0,2), M(3,0)
	double * moments1 = new double[6];
	double * moments2 = new double[6];
	for (int i = 0; i < 6; i++)
	{
		moments1[i] = MAX(moments[i], 1.0e-30);
		moments2[i] = moments[i];
		moments[i] = moments1[i];
	}
		
	
	Double d;
/*	if (-6.0*moments[5]*moments[0]*moments[1]*moments[3]+4.0*moments[5]*pow(moments[1], 3.0)-3.0*pow(moments[3]*moments[1], 2.0)+pow(moments[5]*moments[0], 2.0)+4.0*moments[0]*pow(moments[3], 3.0) == 0.0)
		d = 0.0;
	else*/ if (
36.0*pow(moments[0], 3.0)*pow(moments[3], 5.0)*pow(moments[1], 2.0)
+4.0*pow(moments[0], 5.0)*pow(moments[3], 3.0)*pow(moments[5], 2.0)
-75.0*pow(moments[1]*moments[3], 4.0)*pow(moments[0], 2.0)
+8.0*pow(moments[5]*moments[1], 3.0)*pow(moments[0], 4.0)
+20.0*pow(moments[5]*moments[0], 2.0)*pow(moments[1], 6.0)
+52.0*pow(moments[3], 3.0)*pow(moments[1], 6.0)*moments[0]
+pow(moments[5], 4.0)*pow(moments[0], 6.0)
-12.0*pow(moments[3], 2.0)*pow(moments[1], 8.0)
+16.0*moments[5]*pow(moments[1], 9.0)
-24.0*pow(moments[0]*moments[3], 4.0)*moments[5]*moments[1]
-20.0*pow(moments[0]*moments[3]*moments[1], 3.0)*moments[5]
+42.0*pow(moments[5]*moments[1]*moments[3], 2.0)*pow(moments[0], 4.0)
+96.0*pow(moments[1], 5.0)*pow(moments[0]*moments[3], 2.0)*moments[5]
-12.0*pow(moments[5], 3.0)*pow(moments[0], 5.0)*moments[1]*moments[3]
-60.0*pow(moments[5], 2.0)*pow(moments[0], 3.0)*pow(moments[1], 4.0)*moments[3]
-72.0*moments[5]*moments[0]*pow(moments[1], 7.0)*moments[3]
		< 1.0e-75)
		d = 0.5 * moments[0];
	else
		d = (
4.0*pow(moments[0], 2.0)*pow(moments[3], 3.0)
-3.0*pow(moments[1]*moments[3], 2.0)*moments[0]
-6.0*moments[5]*pow(moments[0], 2.0)*moments[1]*moments[3]
+pow(moments[5], 2.0)*pow(moments[0], 3.0)
+4.0*moments[5]*moments[0]*pow(moments[1], 3.0)-sqrt(MAX(
+36.0*pow(moments[0], 3.0)*pow(moments[3], 5.0)*pow(moments[1], 2.0)
+4.0*pow(moments[0], 5.0)*pow(moments[3], 3.0)*pow(moments[5], 2.0)
-75.0*pow(moments[1]*moments[3], 4.0)*pow(moments[0], 2.0)
+8.0*pow(moments[5]*moments[1], 3.0)*pow(moments[0], 4.0)
+20.0*pow(moments[5]*moments[0], 2.0)*pow(moments[1], 6.0)
+52.0*pow(moments[3], 3.0)*pow(moments[1], 6.0)*moments[0]
+pow(moments[5], 4.0)*pow(moments[0], 6.0)
-12.0*pow(moments[3], 2.0)*pow(moments[1], 8.0)
+16.0*moments[5]*pow(moments[1], 9.0)
-24.0*pow(moments[0]*moments[3], 4.0)*moments[5]*moments[1]
-20.0*pow(moments[0]*moments[3]*moments[1], 3.0)*moments[5]
+42.0*pow(moments[5]*moments[1]*moments[3], 2.0)*pow(moments[0], 4.0)
+96.0*pow(moments[1], 5.0)*pow(moments[0]*moments[3], 2.0)*moments[5]
-12.0*pow(moments[5], 3.0)*pow(moments[0], 5.0)*moments[1]*moments[3]
-60.0*pow(moments[5], 2.0)*pow(moments[0], 3.0)*pow(moments[1], 4.0)*moments[3]
-72.0*moments[5]*moments[0]*pow(moments[1], 7.0)*moments[3], 0.0))) / (2.0*(
-6.0*moments[5]*moments[0]*moments[1]*moments[3]
+4.0*moments[5]*pow(moments[1], 3.0)
-3.0*pow(moments[3]*moments[1], 2.0)
+pow(moments[5]*moments[0], 2.0)
+4.0*moments[0]*pow(moments[3], 3.0)));
	Double f;
	if (d == 0.0)
		f = moments[1] / moments[0];
	else
		f = (moments[2]*d + sqrt(MAX(pow(moments[2]*d, 2.0) + d*moments[4]*pow(moments[0], 2.0) - pow(d, 2.0)*moments[0]*moments[4] - d*moments[0]*pow(moments[2], 2.0), 0.0))) / (d * moments[0]);
	Double e;
	if (d == 0.0)
		e = moments[2] / moments[0];
	else
		e = (moments[1]*d + sqrt(MAX(pow(moments[1]*d, 2.0) + d*moments[3]*pow(moments[0], 2.0) - pow(d, 2.0)*moments[0]*moments[3] - d*moments[0]*pow(moments[1], 2.0), 0.0))) / (d * moments[0]);
	Double c = (moments[2] - d*f) / (moments[0] - d);
	Double b = (moments[1] - d*e) / (moments[0] - d);
	Double a = moments[0] - d;

	Double aexp = a;
	Double bexp = pow(b, Index1);
	Double cexp = pow(c, Index2);
	Double dexp = d;
	Double eexp = pow(e, Index1);
	Double fexp = pow(f, Index2);

	for (int i = 0; i < 6; ++i)
	{
		moments[i] = moments2[i];
	}
	
	delete [] moments1;
	delete [] moments2;

/*if (isinf(aexp * bexp * cexp * dexp * eexp * fexp))
{
	cerr << "Moment(" << Index1 << ", " << Index2 << ") is infinity...\n";
}*/
//if ((aexp * bexp * cexp * dexp * eexp * fexp) > 1)
//{
/*if (moments[1] / moments[0] > 20)
{
	cerr << "Moment(0): " << moments[0] << endl;
	cerr << "Moment(1): " << moments[1] << endl;
	cerr << "Moment(2): " << moments[2] << endl;
	cerr << "Moment(3): " << moments[3] << endl;
	cerr << "Moment(4): " << moments[4] << endl;
	cerr << "Moment(5): " << moments[5] << endl;
	cerr << "a, b, c, d, e, f: " << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << endl;
	cerr << "Moment(" << Index1 << ", " << Index2 << "): " << aexp * bexp * cexp + dexp * eexp * fexp << endl;}
//	cin.get();*/
//}
//	delete temp;
//	return aexp * bexp * cexp * dexp * eexp * fexp * moments[0];
//cerr << "Moment(" << Index1 << ", " << Index2 << "): " << aexp * bexp * cexp * dexp * eexp * fexp << endl;
//cin.get();
	return aexp * bexp * cexp + dexp * eexp * fexp;
	return aexp * bexp * cexp * dexp * eexp * fexp;
}
#endif

//Mueller (8/14/07)
Double TSoot::FractionalMoment4(Double Index1, Double Index2, Double *moments)
{
	//M(0,0), M(0,1), M(1,0), M(0,2), M(1,1), M(2,0), M(0,3), M(1,2), M(2,1), M(3,0)
	Double a = pow(MAX(moments[9], 1.0e-30), 1.0/6.0) * pow(MAX(moments[2], 1.0e-30), 1.0/2.0) * pow(MAX(moments[5], 1.0e-30), -1.0/2.0) * pow(MAX(moments[0], 1.0e-30), -1.0/6.0);
	Double b = pow(MAX(moments[6], 1.0e-30), 1.0/6.0) * pow(MAX(moments[1], 1.0e-30), 1.0/2.0) * pow(MAX(moments[3], 1.0e-30), -1.0/2.0) * pow(MAX(moments[0], 1.0e-30), -1.0/6.0);
	Double c = pow(MAX(moments[5], 1.0e-30), 2.0) * pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[9], 1.0e-30), -1.0/2.0) * pow(MAX(moments[2], 1.0e-30), -5.0/2.0);
	Double d = pow(MAX(moments[3], 1.0e-30), 2.0) * pow(MAX(moments[0], 1.0e-30), 1.0) * pow(MAX(moments[6], 1.0e-30), -1.0/2.0) * pow(MAX(moments[1], 1.0e-30), -5.0/2.0);
	Double e = pow(MAX(moments[9], 1.0e-30), 1.0/3.0) * pow(MAX(moments[2], 1.0e-30), 3.0) * pow(MAX(moments[5], 1.0e-30), -3.0/2.0) * pow(MAX(moments[0], 1.0e-30), -11.0/6.0);
	Double f = pow(MAX(moments[6], 1.0e-30), 1.0/3.0) * pow(MAX(moments[1], 1.0e-30), 3.0) * pow(MAX(moments[3], 1.0e-30), -3.0/2.0) * pow(MAX(moments[0], 1.0e-30), -11.0/6.0);
	Double g = MAX(moments[0], 1.0e-30);
	Double h = pow(MAX(moments[2], 1.0e-30), 1.0) * pow(MAX(moments[1], 1.0e-30), 1.0/2.0) * pow(MAX(moments[8], 1.0e-30), 1.0/2.0) * pow(MAX(moments[0], 1.0e-30), -1.0/2.0) * pow(MAX(moments[5], 1.0e-30), -1.0/2.0) * pow(MAX(moments[4], 1.0e-30), -1.0);
	Double i = pow(MAX(moments[1], 1.0e-30), 1.0) * pow(MAX(moments[2], 1.0e-30), 1.0/2.0) * pow(MAX(moments[7], 1.0e-30), 1.0/2.0) * pow(MAX(moments[0], 1.0e-30), -1.0/2.0) * pow(MAX(moments[3], 1.0e-30), -1.0/2.0) * pow(MAX(moments[4], 1.0e-30), -1.0);
	Double j = pow(MAX(moments[0], 1.0e-30), 2.0) * pow(MAX(moments[5], 1.0e-30), 1.0/2.0) * pow(MAX(moments[3], 1.0e-30), 1.0/2.0) * pow(MAX(moments[4], 1.0e-30), 3.0) * pow(MAX(moments[2], 1.0e-30), -5.0/2.0) * pow(MAX(moments[1], 1.0e-30), -5.0/2.0) * pow(MAX(moments[7], 1.0e-30), -1.0/2.0) * pow(MAX(moments[8], 1.0e-30), -1.0/2.0);

	Double aexp = pow(a, Index1*Index1*Index1);
	Double bexp = pow(b, Index2*Index2*Index2);
	Double cexp = pow(c, Index1*Index1);
	Double dexp = pow(d, Index2*Index2);
	Double eexp = pow(e, Index1);
	Double fexp = pow(f, Index2);
	Double gexp = g;
	Double hexp = pow(h, Index1*Index1*Index2);
	Double iexp = pow(i, Index1*Index2*Index2);
	Double jexp = pow(j, Index1*Index2);

	return aexp * bexp * cexp * dexp * eexp * fexp * gexp * hexp * iexp * jexp;
}

Double TSoot::FractionalMoment5(Double Index1, Double Index2, Double *moments)
{
	//M(0,0), M(1,0), M(0,1), M(0,-1), M(1,-1), M(2,-1), M(-1,-1), M(-1,0), M(-1,1), M(-1,2), M(-2,-1), M(-2,0), M(-2,1), M(-2,2), M(-2,3)
	Double a = pow(MAX(moments[3], 1.0e-30), 1.0/4.0) * pow(MAX(moments[4], 1.0e-30), -1.0/6.0) * pow(MAX(moments[5], 1.0e-30), 1.0/24.0) * pow(MAX(moments[6], 1.0e-30), -1.0/6.0) * pow(MAX(moments[10], 1.0e-30), 1.0/24.0);
	Double b = pow(MAX(moments[0], 1.0e-30), -1.0/2.0) * pow(MAX(moments[1], 1.0e-30), 1.0/6.0) * pow(MAX(moments[3], 1.0e-30), 1.0/2.0) * pow(MAX(moments[4], 1.0e-30), -1.0/3.0) * pow(MAX(moments[5], 1.0e-30), 1.0/12.0) * pow(MAX(moments[6], 1.0e-30), -1.0/3.0) * pow(MAX(moments[7], 1.0e-30), 1.0/2.0) * pow(MAX(moments[10], 1.0e-30), 1.0/12.0) * pow(MAX(moments[11], 1.0e-30), -1.0/6.0);
	Double c = pow(MAX(moments[0], 1.0e-30), -1.0) * pow(MAX(moments[1], 1.0e-30), 1.0/2.0) * pow(MAX(moments[3], 1.0e-30), -1.0/4.0) * pow(MAX(moments[4], 1.0e-30), 1.0/6.0) * pow(MAX(moments[5], 1.0e-30), -1.0/24.0) * pow(MAX(moments[6], 1.0e-30), 1.0/6.0) * pow(MAX(moments[7], 1.0e-30), 1.0/2.0) * pow(MAX(moments[10], 1.0e-30), -1.0/24.0);
	Double d = pow(MAX(moments[0], 1.0e-30), 1.0/2.0) * pow(MAX(moments[1], 1.0e-30), 1.0/3.0) * pow(MAX(moments[3], 1.0e-30), -1.0/2.0) * pow(MAX(moments[4], 1.0e-30), 1.0/3.0) * pow(MAX(moments[5], 1.0e-30), -1.0/12.0) * pow(MAX(moments[6], 1.0e-30), 1.0/3.0) * pow(MAX(moments[7], 1.0e-30), -1.0) * pow(MAX(moments[10], 1.0e-30), -1.0/12.0) * pow(MAX(moments[11], 1.0e-30), 1.0/6.0);
	Double e = pow(MAX(moments[0], 1.0e-30), -1.0/2.0) * pow(MAX(moments[1], 1.0e-30), 1.0/6.0) * pow(MAX(moments[3], 1.0e-30), 1.0/2.0) * pow(MAX(moments[4], 1.0e-30), -1.0/6.0) * pow(MAX(moments[6], 1.0e-30), -1.0/2.0) * pow(MAX(moments[7], 1.0e-30), 1.0/2.0) * pow(MAX(moments[10], 1.0e-30), 1.0/6.0) * pow(MAX(moments[11], 1.0e-30), -1.0/6.0);
	Double f = pow(MAX(moments[0], 1.0e-30), -1.0/2.0) * pow(MAX(moments[2], 1.0e-30), 1.0/4.0) * pow(MAX(moments[3], 1.0e-30), 1.0/4.0) * pow(MAX(moments[6], 1.0e-30), -1.0/2.0) * pow(MAX(moments[7], 1.0e-30), 1.0) * pow(MAX(moments[8], 1.0e-30), -1.0/2.0) * pow(MAX(moments[10], 1.0e-30), 1.0/4.0) * pow(MAX(moments[11], 1.0e-30), -1.0/2.0) * pow(MAX(moments[12], 1.0e-30), 1.0/4.0);
	Double g = pow(MAX(moments[6], 1.0e-30), -1.0/6.0) * pow(MAX(moments[7], 1.0e-30), 1.0/2.0) * pow(MAX(moments[8], 1.0e-30), -1.0/2.0) * pow(MAX(moments[9], 1.0e-30), 1.0/6.0) * pow(MAX(moments[10], 1.0e-30), 1.0/6.0) * pow(MAX(moments[11], 1.0e-30), -1.0/2.0) * pow(MAX(moments[12], 1.0e-30), 1.0/2.0) * pow(MAX(moments[13], 1.0e-30), -1.0/6.0);
	Double h = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[1], 1.0e-30), 1.0/2.0) * pow(MAX(moments[2], 1.0e-30), 1.0/4.0) * pow(MAX(moments[3], 1.0e-30), 5.0/4.0) * pow(MAX(moments[4], 1.0e-30), -1.0/2.0) * pow(MAX(moments[6], 1.0e-30), -1.0) * pow(MAX(moments[7], 1.0e-30), 3.0/2.0) * pow(MAX(moments[8], 1.0e-30), -1.0/2.0) * pow(MAX(moments[10], 1.0e-30), 1.0/4.0) * pow(MAX(moments[11], 1.0e-30), -1.0/2.0) * pow(MAX(moments[12], 1.0e-30), 1.0/4.0);
	Double i = pow(MAX(moments[0], 1.0e-30), -1.0) * pow(MAX(moments[1], 1.0e-30), 1.0/3.0) * pow(MAX(moments[2], 1.0e-30), 3.0/4.0) * pow(MAX(moments[3], 1.0e-30), 1.0/4.0) * pow(MAX(moments[4], 1.0e-30), -1.0/3.0) * pow(MAX(moments[6], 1.0e-30), 1.0/6.0) * pow(MAX(moments[7], 1.0e-30), 1.0/2.0) * pow(MAX(moments[8], 1.0e-30), -1.0/2.0) * pow(MAX(moments[9], 1.0e-30), -1.0/6.0) * pow(MAX(moments[10], 1.0e-30), -1.0/12.0) * pow(MAX(moments[11], 1.0e-30), 1.0/6.0) * pow(MAX(moments[12], 1.0e-30), -1.0/4.0) * pow(MAX(moments[13], 1.0e-30), 1.0/6.0);
	Double j = pow(MAX(moments[0], 1.0e-30), -3.0/2.0) * pow(MAX(moments[2], 1.0e-30), 3.0/4.0) * pow(MAX(moments[3], 1.0e-30), 3.0/4.0) * pow(MAX(moments[6], 1.0e-30), -1.0) * pow(MAX(moments[7], 1.0e-30), 2.0) * pow(MAX(moments[8], 1.0e-30), -1.0) * pow(MAX(moments[10], 1.0e-30), 1.0/4.0) * pow(MAX(moments[11], 1.0e-30), -1.0/2.0) * pow(MAX(moments[12], 1.0e-30), 1.0/4.0);
	Double k = pow(MAX(moments[2], 1.0e-30), 1.0/2.0) * pow(MAX(moments[3], 1.0e-30), -1.0/2.0) * pow(MAX(moments[6], 1.0e-30), 1.0/3.0) * pow(MAX(moments[7], 1.0e-30), -1.0) * pow(MAX(moments[8], 1.0e-30), 1.0) * pow(MAX(moments[9], 1.0e-30), -1.0/3.0) * pow(MAX(moments[10], 1.0e-30), -1.0/12.0) * pow(MAX(moments[11], 1.0e-30), 1.0/6.0) * pow(MAX(moments[13], 1.0e-30), -1.0/6.0) * pow(MAX(moments[14], 1.0e-30), 1.0/12.0);
	Double l = pow(MAX(moments[0], 1.0e-30), -1.0) * pow(MAX(moments[2], 1.0e-30), 1.0/2.0) * pow(MAX(moments[3], 1.0e-30), 1.0/2.0) * pow(MAX(moments[10], 1.0e-30), -1.0/24.0) * pow(MAX(moments[11], 1.0e-30), 1.0/6.0) * pow(MAX(moments[12], 1.0e-30), -1.0/4.0) * pow(MAX(moments[13], 1.0e-30), 1.0/6.0) * pow(MAX(moments[14], 1.0e-30), -1.0/24.0);
	Double m = pow(MAX(moments[6], 1.0e-30), -1.0/3.0) * pow(MAX(moments[7], 1.0e-30), 1.0) * pow(MAX(moments[8], 1.0e-30), -1.0) * pow(MAX(moments[9], 1.0e-30), 1.0/3.0) * pow(MAX(moments[10], 1.0e-30), 1.0/12.0) * pow(MAX(moments[11], 1.0e-30), -1.0/6.0) * pow(MAX(moments[13], 1.0e-30), 1.0/6.0) * pow(MAX(moments[14], 1.0e-30), -1.0/12.0);
	Double n = pow(MAX(moments[10], 1.0e-30), 1.0/24.0) * pow(MAX(moments[11], 1.0e-30), -1.0/6.0) * pow(MAX(moments[12], 1.0e-30), 1.0/4.0) * pow(MAX(moments[13], 1.0e-30), -1.0/6.0) * pow(MAX(moments[14], 1.0e-30), 1.0/24.0);
	Double q = MAX(moments[0], 1.0e-30);

	Double aexp = pow(a, Index1*Index1*Index1*Index1/20.0);
	Double bexp = pow(b, Index1*Index1*Index1/20.0);
	Double cexp = pow(c, Index1*Index1/20.0);
	Double dexp = pow(d, Index1/20.0);
	Double eexp = pow(e, Index1*Index1*Index1*Index2/20.0);
	Double fexp = pow(f, Index1*Index1*Index2*Index2/20.0);
	Double gexp = pow(g, Index1*Index2*Index2*Index2/20.0);
	Double hexp = pow(h, Index1*Index1*Index2/20.0);
	Double iexp = pow(i, Index1*Index2/20.0);
	Double jexp = pow(j, Index1*Index2*Index2/20.0);
	Double kexp = pow(k, Index2/20.0);
	Double lexp = pow(l, Index2*Index2/20.0);
	Double mexp = pow(m, Index2*Index2*Index2/20.0);
	Double nexp = pow(n, Index2*Index2*Index2*Index2/20.0);
	Double qexp = pow(q, 1.0/20.0);

/*if (Index1==-2.5 && Index2==3)
{
	cerr << "aexp: " << aexp << endl;
	cerr << "bexp: " << bexp << endl;
	cerr << "cexp: " << cexp << endl;
	cerr << "dexp: " << dexp << endl;
	cerr << "eexp: " << eexp << endl;
	cerr << "fexp: " << fexp << endl;
	cerr << "gexp: " << gexp << endl;
	cerr << "hexp: " << hexp << endl;
	cerr << "iexp: " << iexp << endl;
	cerr << "jexp: " << jexp << endl;
	cerr << "kexp: " << kexp << endl;
	cerr << "lexp: " << lexp << endl;
	cerr << "mexp: " << mexp << endl;
	cerr << "nexp: " << nexp << endl;
	cerr << "qexp: " << qexp << endl;
	cerr << "Product: " << aexp * bexp * cexp * dexp * eexp * fexp * gexp * hexp * iexp * jexp * kexp * lexp * mexp * nexp * qexp * aexp * bexp * cexp * dexp * eexp * fexp * gexp * hexp * iexp * jexp * kexp * lexp * mexp * nexp * qexp<< endl;
	cin.get();
}*/

	return pow(aexp * bexp * cexp * dexp * eexp * fexp * gexp * hexp * iexp * jexp * kexp * lexp * mexp * nexp * qexp, 20.0);
}

Double TSoot::FractionalMoment2( Double fracIndex, int first, int second, Double *moments )
{
/*	if ( fracIndex < 0.0 || fracIndex > 1.0 ) {*/
/*		return 0.0;*/
/*	}*/
	if ( moments[first] < SMALLSOOT || moments[second] < SMALLSOOT ) {
		return 0.0;
	}

#ifdef CUTFRACMOM
#	ifdef CUTLOWER
	return MAX( exp( log( MAX( moments[first], SMALLSOOT ) )
						* GetAlphaI2( first, fracIndex, first, second )
				+ log( MAX( moments[second], SMALLSOOT ) ) 
						* GetAlphaI2( second, fracIndex, first, second ) )
			, MAX( moments[first], SMALLSOOT ) );
#	else
	return MAX( MIN( exp( log( MAX( moments[first], SMALLSOOT ) )
						* GetAlphaI2( first, fracIndex, first, second )
				+ log( MAX( moments[second], SMALLSOOT ) ) 
						* GetAlphaI2( second, fracIndex, first, second ) )
			, moments[second] ), MAX( moments[first], SMALLSOOT ) );
#	endif
#else
	return exp( log( MAX( moments[first], SMALLSOOT ) ) 
						* GetAlphaI2( first, fracIndex, first, second )
					+ log( MAX( moments[second], SMALLSOOT ) ) 
						* GetAlphaI2( second, fracIndex, first, second ) );
#endif
}

//Mueller (7/20/07)
Double TSoot::FracMom2(Double Index1, Double Index2, Double *moments)
{
	if (fNSootMoments==3)
	{
		return FractionalMoment2(Index1, Index2, moments);
	}
	else if (fNSootMoments==6)
	{
		return FractionalMoment3(Index1, Index2, moments);
	}
	else if (fNSootMoments==10)
	{
		return FractionalMoment4(Index1, Index2, moments);
	}
	else if (fNSootMoments==15)
	{
		return FractionalMoment5(Index1, Index2, moments);
	}
	else
	{
		return 0.0;
	}
	return 0.0;
}

Double TSoot::FracMom2( Double fracIndex, Double *moments )
{
#ifdef NEWINTERPOLATION
	int	second;

	if ( fracIndex < 0.0 ) {
		second = 1;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		second = fNSootMoments-1;
	}
	else {
		second = ( int )fracIndex+1;
	}

	return FractionalMoment2( fracIndex, 0, second, moments );
#else
	int	first;
	
	if ( fracIndex < 0.0 ) {
		first = 0;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		first = fNSootMoments-2;
	}
	else {
		first = ( int )fracIndex;
	}

	return FractionalMoment2( fracIndex, first, first+1, moments );
#endif
}

/*Double TSoot::GetAlphaI( int i, Double fracIndex )
{
	const Double oneOverSix = -1.0 / 6.0;
	switch ( i ) {
		case 0:
			return oneOverSix * ( 1.0 - fracIndex ) * ( fracIndex - 2.0 ) * ( fracIndex - 3.0 );
		case 1:
			return 0.5 * fracIndex * ( fracIndex - 2.0 ) * ( fracIndex - 3.0 );
		case 2:
			return 0.5 * fracIndex * ( 1.0 - fracIndex ) * ( fracIndex - 3.0 );
		case 3:
			return oneOverSix * fracIndex * ( fracIndex - 1.0 ) * ( fracIndex - 2.0 );
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
}
*/

//Mueller (7/19/07)
Double TSoot::GetAlphaI2(int i, Double fracIndex, int first, int second)
{
	if (i == first)
	{
		return (fracIndex-second) / (first-second);
	}
	else if (i == second)
	{
		return (fracIndex-first) / (second-first);
	}
	else
	{
		return 0.0;
	}
}

//Mueller (8/7/07)
Double TSoot::GetAlphaI2(int i, Double fracIndex, int first, int second, int third)
{
	if (i==first)
	{
		return ((fracIndex-second) / (first-second)) * ((fracIndex-third) / (first-third)) ;
	}
	else if (i==second)
	{
		return ((fracIndex-first) / (second-first)) * ((fracIndex-third) / (second-third));
	}
	else if (i==third)
	{
		return ((fracIndex-first) / (third-first)) * ((fracIndex-second) / (third-second));
	}
	else
	{
		return 0.0;
	}
}

Double TSoot::GetAlphaI2( int i, Double fracIndex )
{	
#ifdef NEWINTERPOLATION
	int	second;

	if ( fracIndex < 0.0 ) {
		second = 1;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		second = fNSootMoments-1;
	}
	else {
		second = ( int )fracIndex+1;
	}

	return GetAlphaI2( i, fracIndex, 0, second );
#else
	int	first;

	if ( fracIndex < 0.0 ) {
		first = 0;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		first = fNSootMoments-2;
	}
	else {
		first = ( int )fracIndex;
	}

	return GetAlphaI2( i, fracIndex, first, first+1 );
#endif
}

void TSoot::ComputeFractionalMomentsPAH( Double *moments )
{
	Double	*fracMom = fFracMomPAH->vec;

	fracMom[km1_2] = FracMom2( -0.5, moments );
	fracMom[k01_2] = FracMom2( 0.5, moments );
	fracMom[k03_2] = FracMom2( 1.5, moments );
	fracMom[k05_2] = FracMom2( 2.5, moments );

	fracMom[k19_6] = FracMom2( 19.0 / 6.0, moments );
	fracMom[k17_6] = FracMom2( 17.0 / 6.0, moments );
	fracMom[k13_6] = FracMom2( 13.0 / 6.0, moments );
	fracMom[k11_6] = FracMom2( 11.0 / 6.0, moments );

	fracMom[k07_6] = FracMom2( 7.0 / 6.0, moments );
	fracMom[k05_6] = FracMom2( 5.0 / 6.0, moments );
	fracMom[k01_6] = FracMom2( 1.0 / 6.0, moments );
	fracMom[km1_6] = FracMom2( -1.0 / 6.0, moments );
}

void TSoot::ComputeFractionalMoments( Double *moments )
{
	Double	*fracMom = fFracMom->vec;
	
	fracMom[km1_2] = FracMom2( -0.5, moments );
	fracMom[k01_2] = FracMom2( 0.5, moments );
	fracMom[k03_2] = FracMom2( 1.5, moments );
	fracMom[k05_2] = FracMom2( 2.5, moments );

	fracMom[k19_6] = FracMom2( 19.0 / 6.0, moments );
	fracMom[k17_6] = FracMom2( 17.0 / 6.0, moments );
	fracMom[k13_6] = FracMom2( 13.0 / 6.0, moments );
	fracMom[k11_6] = FracMom2( 11.0 / 6.0, moments );

	fracMom[k07_6] = FracMom2( 7.0 / 6.0, moments );
	fracMom[k05_6] = FracMom2( 5.0 / 6.0, moments );
	fracMom[k01_6] = FracMom2( 1.0 / 6.0, moments );
	fracMom[km1_6] = FracMom2( -1.0 / 6.0, moments );
	
	fracMom[km1_24] = FracMom2( -1.0 / 24.0, moments );
	fracMom[k23_24] = FracMom2( 23.0 / 24.0, moments );
	fracMom[k47_24] = FracMom2( 47.0 / 24.0, moments );
	fracMom[k71_24] = FracMom2( 71.0 / 24.0, moments );

	fracMom[km1_3] = FracMom2( -1.0 / 3.0, moments );
	fracMom[k02_3] = FracMom2( 2.0 / 3.0, moments );
	fracMom[k05_3] = FracMom2( 5.0 / 3.0, moments );
	fracMom[k08_3] = FracMom2( 8.0 / 3.0, moments );

}

Double TSoot::GetDerivFracMom( int l, Double which, Double *moments )
{
	// returns dM_which/dM_l
	return GetAlphaI2( l, which ) * FracMom2( which, moments ) / moments[l];
}

Double TSoot::GetDerivFracMom( int l, FracMoments nMom, Double *moments )
{
	// returns dM_r/dM_l

	Double	*fm = fFracMom->vec;

	switch( nMom ) {
		case km1_2:
			return GetAlphaI2( l, -1.0/2.0 ) * fm[nMom] / moments[l];
		case k01_2:
			return GetAlphaI2( l, 1.0/2.0 ) * fm[nMom] / moments[l];
		case k05_6:
			return GetAlphaI2( l, 5.0/6.0 ) * fm[nMom] / moments[l];
		case k01_6:
			return GetAlphaI2( l, 1.0/6.0 ) * fm[nMom] / moments[l];
		case km1_6:
			return GetAlphaI2( l, -1.0/6.0 ) * fm[nMom] / moments[l];

		case k03_2:
			return GetAlphaI2( l, 1.5 ) * fm[nMom] / moments[l];
		case k11_6:
			return GetAlphaI2( l, 11.0/6.0 ) * fm[nMom] / moments[l];
		case k07_6:
			return GetAlphaI2( l, 7.0/6.0 ) * fm[nMom] / moments[l];

		case k05_2:
			return GetAlphaI2( l, 2.5 ) * fm[nMom] / moments[l];
		case k19_6:
			return GetAlphaI2( l, 19.0/6.0 ) * fm[nMom] / moments[l];
		case k17_6:
			return GetAlphaI2( l, 17.0/6.0 ) * fm[nMom] / moments[l];
		case k13_6:
			return GetAlphaI2( l, 13.0/6.0) * fm[nMom] / moments[l];

		case km1_24:
			return GetAlphaI2( l, -1.0/24.0 ) * fm[nMom] / moments[l];
		case k23_24:
			return GetAlphaI2( l, 23.0/24.0 ) * fm[nMom] / moments[l];
		case k47_24:
			return GetAlphaI2( l, 47.0/24.0 ) * fm[nMom] / moments[l];
		case k71_24:
			return GetAlphaI2( l, 71.0/24.0 ) * fm[nMom] / moments[l];

		case km1_3:
			return GetAlphaI2( l, -1.0/3.0 ) * fm[nMom] / moments[l];
		case k02_3:
			return GetAlphaI2( l, 2.0/3.0 ) * fm[nMom] / moments[l];
		case k05_3:
			return GetAlphaI2( l, 5.0/3.0 ) * fm[nMom] / moments[l];
		case k08_3:
			return GetAlphaI2( l, 8.0/3.0 ) * fm[nMom] / moments[l];

		default:
			cerr << "#error: no need to compute derivative of alpha_" << (int)nMom << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::GetDerivPhi( int l, Phi which, Double *moments )
{
	Double	*phi = fPhi->vec;
	Double	*fm = fFracMom->vec;

	switch( which ) {
		case kh00:
			return 0.5 * phi[kh00] 
					* ( GetDerivPhi( l, k000, moments ) / phi[k000]
						+ GetDerivPhi( l, k100, moments ) / phi[k100] );
		case kh11:
			return 0.5 * phi[kh11] 
					* ( GetDerivPhi( l, k011, moments ) / phi[k011]
						+ GetDerivPhi( l, k111, moments ) / phi[k111] );
		case kh12:
			return 0.5 * phi[kh12] 
					* ( GetDerivPhi( l, k012, moments ) / phi[k012]
						+ GetDerivPhi( l, k112, moments ) / phi[k112] );
		case k000:
			return    2.0 * ( fm[k01_6] * GetDerivFracMom( l, km1_2, moments ) 
							+ fm[km1_2] * GetDerivFracMom( l, k01_6, moments )
							+ 2.0 * fm[km1_6] * GetDerivFracMom( l, km1_6, moments ) );
		case k100:
			return    2.0 * ( fm[k07_6] * GetDerivFracMom( l, km1_2, moments ) 
							+ fm[km1_2] * GetDerivFracMom( l, k07_6, moments )
					+ 2.0 * ( fm[k05_6] * GetDerivFracMom( l, km1_6, moments )
							+ fm[km1_6] * GetDerivFracMom( l, k05_6, moments ) )
							+ fm[k01_2] * GetDerivFracMom( l, k01_6, moments )
							+ fm[k01_6] * GetDerivFracMom( l, k01_2, moments ) );
		case k011:
			return    2.0 * ( fm[k07_6] * GetDerivFracMom( l, k01_2, moments ) 
							+ fm[k01_2] * GetDerivFracMom( l, k07_6, moments )
					  + 2.0 * fm[k05_6] * GetDerivFracMom( l, k05_6, moments ) );
		case k111:
			return    2.0 * ( fm[k13_6] * GetDerivFracMom( l, k01_2, moments ) 
							+ fm[k01_2] * GetDerivFracMom( l, k13_6, moments )
					+ 2.0 * ( fm[k11_6] * GetDerivFracMom( l, k05_6, moments )
							+ fm[k05_6] * GetDerivFracMom( l, k11_6, moments ) )
							+ fm[k03_2] * GetDerivFracMom( l, k07_6, moments )
							+ fm[k07_6] * GetDerivFracMom( l, k03_2, moments ) );
		case k012:
			return            fm[k07_6] * GetDerivFracMom( l, k03_2, moments ) 
							+ fm[k03_2] * GetDerivFracMom( l, k07_6, moments )
					+ 2.0 * ( fm[k05_6] * GetDerivFracMom( l, k11_6, moments )
							+ fm[k11_6] * GetDerivFracMom( l, k05_6, moments ) )
							+ fm[k01_2] * GetDerivFracMom( l, k13_6, moments )
							+ fm[k13_6] * GetDerivFracMom( l, k01_2, moments );
		case k112:
			return 		      fm[k13_6] * GetDerivFracMom( l, k03_2, moments ) 
							+ fm[k03_2] * GetDerivFracMom( l, k13_6, moments )
			  + 2.0 * ( 2.0 * fm[k11_6] * GetDerivFracMom( l, k11_6, moments )
							+ fm[k05_6] * GetDerivFracMom( l, k17_6, moments )
							+ fm[k17_6] * GetDerivFracMom( l, k05_6, moments ) )

							+ fm[k03_2] * GetDerivFracMom( l, k13_6, moments )
							+ fm[k13_6] * GetDerivFracMom( l, k03_2, moments )
							+ fm[k07_6] * GetDerivFracMom( l, k05_2, moments )
							+ fm[k05_2] * GetDerivFracMom( l, k07_6, moments )
							+ fm[k01_2] * GetDerivFracMom( l, k19_6, moments )
							+ fm[k19_6] * GetDerivFracMom( l, k01_2, moments );
		default:
			cerr << "#error: no need to compute derivative of phi_" << (int)which << NEWL;
			exit( 2 );
	}
	return 0.0;
}

void TSoot::ComputePhi( Double temp )
{
	Double	c = GetC( temp );
	Double	*fm = fFracMom->vec;
	Double	*phi = fPhi->vec;

	phi[k000] = ( fm[k01_6] * fm[km1_2] + fm[km1_6] * fm[km1_6] ) * 2.0;
	phi[k100] = ( fm[k07_6] * fm[km1_2] + fm[k05_6] * fm[km1_6] * 2.0 
									+ fm[k01_2] * fm[k01_6] ) * 2.0;
	
	phi[k011] = ( fm[k07_6] * fm[k01_2] + fm[k05_6] * fm[k05_6] ) * 2.0;
	phi[k111] = ( fm[k13_6] * fm[k01_2] + fm[k11_6] * fm[k05_6] * 2.0 
									+ fm[k03_2] * fm[k07_6] ) * 2.0;
	
	phi[k012] =   fm[k07_6] * fm[k03_2] + fm[k05_6] * fm[k11_6] * 2.0 
									+ fm[k01_2] * fm[k13_6];
	phi[k112] =   fm[k13_6] * fm[k03_2] 
					+ ( fm[k11_6] * fm[k11_6] + fm[k05_6] * fm[k17_6] ) * 2.0
				 	+ fm[k03_2] * fm[k13_6] + fm[k07_6] * fm[k05_2] + fm[k01_2] * fm[k19_6];
	
	phi[kh00] = sqrt( phi[k000] * phi[k100] ) * c;
	phi[kh11] = sqrt( phi[k011] * phi[k111] ) * c;
	phi[kh12] = sqrt( phi[k012] * phi[k112] ) * c;
	
/*	cerr << "phi[kh00] = " << phi[kh00] << TAB
		 << "phi[kh11] = " << phi[kh11] << TAB
		 << "phi[kh12] = " << phi[kh12] << NEWL;*/
}

/*Double TSoot::GetPhi(Double w, Double y, int x, int z, Double temp, Double Df, Double *moments1, Double *moments2, Double exp)
{
	Double	phi0, phi1, phi2, phi3;
	Double	c = GetC(temp, Df);

	Double av = 1.0 - (2.0 / Df);
	Double as = (3.0 / Df) - 1.0;
	
	phi0 = FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) * FracMom2(x-0.5, z, moments2) + 2.0 * FracMom2(av+w-0.5, as+y, moments1) * FracMom2(av+x-0.5, as+z, moments2) + FracMom2(w-0.5, y, moments1) * FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2);

	phi1 = FracMom2(2.0*av+w+0.5, 2.0*as+y, moments1) * FracMom2(x-0.5, z, moments2) + 2.0 * FracMom2(av+w+0.5, as+y, moments1) * FracMom2(av+x-0.5, as+z, moments2) + FracMom2(w+0.5, y, moments1) * FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2) + FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) * FracMom2(x+0.5, z, moments2) + 2.0 * FracMom2(av+w-0.5, as+y, moments1) * FracMom2(av+x+0.5, as+z, moments2) + FracMom2(w-0.5, y, moments1) * FracMom2(2.0*av+x+0.5, 2.0*as+z, moments2);

	phi2 = FracMom2(2.0*av+w+1.5, 2.0*as+y, moments1) * FracMom2(x-0.5, z, moments2) + 2.0 * FracMom2(av+w+1.5, as+y, moments1) * FracMom2(av+x-0.5, as+z, moments2) + FracMom2(w+1.5, y, moments1) * FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2) + 2.0 * FracMom2(2.0*av+w+0.5, 2.0*as+y, moments1) * FracMom2(x+0.5, z, moments2) + 4.0 * FracMom2(av+w+0.5, as+y, moments1) * FracMom2(av+x+0.5, as+z, moments2) + 2.0 * FracMom2(w+0.5, y, moments1) * FracMom2(2.0*av+x+0.5, 2.0*as+z, moments2) + FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) * FracMom2(x+1.5, z, moments2) + 2.0 * FracMom2(av+w-0.5, as+y, moments1) * FracMom2(av+x+1.5, as+z, moments2) + FracMom2(w-0.5, y, moments1) * FracMom2(2.0*av+x+1.5, 2.0*as+z, moments2);

	phi3 = FracMom2(2.0*av+w+2.5, 2.0*as+y, moments1) * FracMom2(x-0.5, z, moments2) + 2.0 * FracMom2(av+w+2.5, as+y, moments1) * FracMom2(av+x-0.5, as+z, moments2) + FracMom2(w+2.5, y, moments1) * FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2) + 3.0 * FracMom2(2.0*av+w+1.5, 2.0*as+y, moments1) * FracMom2(x+0.5, z, moments2) + 6.0 * FracMom2(av+w+1.5, as+y, moments1) * FracMom2(av+x+0.5, as+z, moments2) + 3.0 * FracMom2(w+1.5, y, moments1) * FracMom2(2.0*av+x+0.5, 2.0*as+z, moments2) + 3.0 * FracMom2(2.0*av+w+0.5, 2.0*as+y, moments1) * FracMom2(x+1.5, z, moments2) + 6.0 * FracMom2(av+w+0.5, as+y, moments1) * FracMom2(av+x+1.5, as+z, moments2) + 3.0 * FracMom2(w+0.5, y, moments1) * FracMom2(2.0*av+x+1.5, 2.0*as+z, moments2) + FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) * FracMom2(x+2.5, z, moments2) + 2.0 * FracMom2(av+w-0.5, as+y, moments1) * FracMom2(av+x+2.5, as+z, moments2) + FracMom2(w-0.5, y, moments1) * FracMom2(2.0*av+x+2.5, 2.0*as+z, moments2);
if (isnan(phi3) || isnan(phi2) || isnan(phi1) || isnan(phi0))
{
	cerr << "One of the phi's is NAN...\n";
	cin.get();
}
	cerr << "Phi0: " << phi0 << endl;
	cerr << "Phi1: " << phi1 << endl;
	cerr << "Phi2: " << phi2 << endl;
	cerr << "Phi3: " << phi3 << endl;

	if (fNSootMoments==3)
	{
		return pow(phi0, 1.0-exp) * pow(phi1, exp) * c;
	}
	else if (fNSootMoments==6)
	{
		return pow(phi0, (1.0-exp)*(1.0-0.5*exp)) * pow(phi1, exp*(2.0-exp)) * pow(phi2, 0.5*exp*(exp-1.0)) * c;
	}
	else if (fNSootMoments==10)
	{
		return pow(phi0, (1.0-exp)*(1.0-0.5*exp)*(1.0-(1.0/3.0)*exp)) * pow(phi1, exp*(2.0-exp)*(1.5-0.5*exp)) * pow(phi2, 0.5*exp*(exp-1.0)*(3.0-exp)) * pow(phi3, (1.0/3.0)*exp*(0.5*exp-0.5)*(exp-2.0)) * c;
	}
}*/

Double TSoot::GetPhi(Double temp, Double *moments1, Double *moments2, Double exp)
{
	Double	phi0, phi1, phi2, phi3;
	Double	c = GetC(temp, 3.0);

/*if (moments1[0] < 0 || moments1[1] < 0 || moments1[2] < 0 || moments1[3] < 0 || moments1[4] < 0 ||moments1[5] < 0)
{
	cerr << "Volume: " << moments1[2] << endl;
	cerr << "Volume 2/3: " << moments1[4] << endl;
	cerr << "M00: " << moments1[0] << endl;
	cerr << " 7.0/6.0: " << FracMom2( 7.0/6.0, 0, moments1) << endl;
	cerr << " 5.0/6.0: " << FracMom2( 5.0/6.0, 0, moments1) << endl;
	cerr << " 1.0/2.0: " << FracMom2( 1.0/2.0, 0, moments1) << endl;
	cerr << " 1.0/6.0: " << FracMom2( 1.0/6.0, 0, moments1) << endl;
	cerr << "-1.0/6.0: " << FracMom2(-1.0/6.0, 0, moments1) << endl;
	cerr << "-1.0/2.0: " << FracMom2(-1.0/2.0, 0, moments1) << endl;
	cin.get();
}*/
	
	phi0 = FracMom2(1.0/6.0, 0, moments1) * FracMom2(-0.5, 0, moments2) + 2.0 * FracMom2(-1.0/6.0, 0, moments1) * FracMom2(-1.0/6.0, 0, moments2) + FracMom2(-0.5, 0, moments1) * FracMom2(1.0/6.0, 0, moments2);

	phi1 = FracMom2(7.0/6.0, 0, moments1) * FracMom2(-0.5, 0, moments2) + 2.0 * FracMom2(5.0/6.0, 0, moments1) * FracMom2(-1.0/6.0, 0, moments2) + FracMom2(0.5, 0, moments1) * FracMom2(1.0/6.0, 0, moments2) + FracMom2(1.0/6.0, 0, moments1) * FracMom2(0.5, 0, moments2) + 2.0 * FracMom2(-1.0/6.0, 0, moments1) * FracMom2(5.0/6.0, 0, moments2) + FracMom2(-0.5, 0, moments1) * FracMom2(7.0/6.0, 0, moments2);

/*	if (fNSootMoments > 3)
	phi2 = FracMom2(13.0/6.0, 0, moments1) * FracMom2(-0.5, 0, moments2) + 2.0 * FracMom2(11.0/6.0, 0, moments1) * FracMom2(-1.0/6.0, 0, moments2) + FracMom2(1.5, 0, moments1) * FracMom2(1.0/6.0, 0, moments2) + 2.0 * FracMom2(7.0/6.0, 0, moments1) * FracMom2(0.5, 0, moments2) + 4.0 * FracMom2(5.0/6.0, 0, moments1) * FracMom2(5.0/6.0, 0, moments2) + 2.0 * FracMom2(0.5, 0, moments1) * FracMom2(7.0/6.0, 0, moments2) + FracMom2(1.0/6.0, 0, moments1) * FracMom2(1.5, 0, moments2) + 2.0 * FracMom2(-1.0/6.0, 0, moments1) * FracMom2(11.0/6.0, 0, moments2) + FracMom2(-0.5, 0, moments1) * FracMom2(13.0/6.0, 0, moments2);

	if (fNSootMoments > 6)
	phi3 = FracMom2(19.0/6.0, 0, moments1) * FracMom2(-0.5, 0, moments2) + 2.0 * FracMom2(17.0/6.0, 0, moments1) * FracMom2(-1.0/6.0, 0, moments2) + FracMom2(2.5, 0, moments1) * FracMom2(1.0/6.0, 0, moments2) + 3.0 * FracMom2(13.0/6.0, 0, moments1) * FracMom2(0.5, 0, moments2) + 6.0 * FracMom2(11.0/6.0, 0, moments1) * FracMom2(5.0/6.0, 0, moments2) + 3.0 * FracMom2(1.5, 0, moments1) * FracMom2(7.0/6.0, 0, moments2) + 3.0 * FracMom2(7.0/6.0, 0, moments1) * FracMom2(1.5, 0, moments2) + 6.0 * FracMom2(5.0/6.0, 0, moments1) * FracMom2(11.0/6.0, 0, moments2) + 3.0 * FracMom2(0.5, 0, moments1) * FracMom2(13.0/6.0, 0, moments2) + FracMom2(1.0/6.0, 0, moments1) * FracMom2(2.5, 0, moments2) + 2.0 * FracMom2(-1.0/6.0, 0, moments1) * FracMom2(17.0/6.0, 0, moments2) + FracMom2(-0.5, 0, moments1) * FracMom2(19.0/6.0, 0, moments2);*/

/*if (phi0 > 1.0e0 || phi1 > 1.0e0)//|| phi2 > 1.0e0)
{
cerr << "Nucleation Phi's:\n";
cerr << "Phi0: " << phi0 << endl;
cerr << "Phi1: " << phi1 << endl;
//cerr << "Phi2: " << phi2 << endl;
}*/

//	if (fNSootMoments==3)
//	{
		return pow(phi0, 1.0-exp) * pow(phi1, exp) * c;
//	}
/*	else if (fNSootMoments==6)
	{
		return pow(phi0, (1.0-exp)*(1.0-0.5*exp)) * pow(phi1, exp*(2.0-exp)) * pow(phi2, 0.5*exp*(exp-1.0)) * c;
	}
	else if (fNSootMoments==10)
	{
		return pow(phi0, (1.0-exp)*(1.0-0.5*exp)*(1.0-(1.0/3.0)*exp)) * pow(phi1, exp*(2.0-exp)*(1.5-0.5*exp)) * pow(phi2, 0.5*exp*(exp-1.0)*(3.0-exp)) * pow(phi3, (1.0/3.0)*exp*(0.5*exp-0.5)*(exp-2.0)) * c;
	}*/
}

//Mueller (7/19/07)
Double TSoot::GetPhi(Double w, Double y, int x, int z, Double temp, Double Df, Double *moments1, Double *moments2)
{
	Double	phi0, phi1;
	Double	c = GetC(temp, Df);

	Double av = 1.0 - (2.0 / Df);
	Double as = (3.0 / Df) - 1.0;

/*if (moments1[0] < 0 || moments1[1] < 0 || moments1[2] < 0 || moments1[3] < 0 || moments1[4] < 0 ||moments1[5] < 0)
{
	cerr << "Volume: " << moments1[2] << endl;
	cerr << "Volume 2/3: " << moments1[4] << endl;
	cerr << "M00: " << moments1[0] << endl;
	cerr << " 7.0/6.0: " << FracMom2( 7.0/6.0, 0, moments1) << endl;
	cerr << " 5.0/6.0: " << FracMom2( 5.0/6.0, 0, moments1) << endl;
	cerr << " 1.0/2.0: " << FracMom2( 1.0/2.0, 0, moments1) << endl;
	cerr << " 1.0/6.0: " << FracMom2( 1.0/6.0, 0, moments1) << endl;
	cerr << "-1.0/6.0: " << FracMom2(-1.0/6.0, 0, moments1) << endl;
	cerr << "-1.0/2.0: " << FracMom2(-1.0/2.0, 0, moments1) << endl;
	cin.get();
}*/
	
	phi0 =	FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) * FracMom2(x-0.5, z, moments2) + 2.0 * FracMom2(av+w-0.5, as+y, moments1) * FracMom2(av+x-0.5, as+z, moments2) + FracMom2(w-0.5, y, moments1) * FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2);

	phi1 =	FracMom2(2.0*av+w+0.5, 2.0*as+y, moments1) * FracMom2(x-0.5, z, moments2) + 2.0 * FracMom2(av+w+0.5, as+y, moments1) * FracMom2(av+x-0.5, as+z, moments2) + FracMom2(w+0.5, y, moments1) * FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2) + FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) * FracMom2(x+0.5, z, moments2) + 2.0 * FracMom2(av+w-0.5, as+y, moments1) * FracMom2(av+x+0.5, as+z, moments2) + FracMom2(w-0.5, y, moments1) * FracMom2(2.0*av+x+0.5, 2.0*as+z, moments2);

/*if (w == -2 && y == 3)
{
	cerr << "Moment(-2,3): " << moments1[14] << endl;
	cerr << "Phi0: " << phi0 << endl;
	cerr << "Moment(-11/6,3): " << FracMom2(2.0*av+w-0.5, 2.0*as+y, moments1) << endl;
	cerr << "Moment(-1/2,0): " << FracMom2(x-0.5, z, moments2) << endl;
	cerr << "Moment(-13/6,3): " << FracMom2(av+w-0.5, as+y, moments1) << endl;
	cerr << "Moment(-1/6,0): " << FracMom2(av+x-0.5, as+z, moments2) << endl;
	cerr << "Moment(-5/2,3): " << FracMom2(w-0.5, y, moments1) << endl;
	cerr << "Moment(1/6,0): " << FracMom2(2.0*av+x-0.5, 2.0*as+z, moments2) << endl;
	cerr << "Phi1: " << phi1 << endl;
}*/

/*if (phi0 > 1.0e0 || phi1 > 1.0e0)//|| phi2 > 1.0e0)
{
cerr << "Coagulation Phi's:\n";
cerr << "Phi0: " << phi0 << endl;
cerr << "Phi1: " << phi1 << endl;
//cerr << "Phi2: " << phi2 << endl;
}
if (phi1 > 1943 && phi1 < 1944)
{
	cerr << "Volume: " << moments1[2] << endl;
	cerr << "Volume 2/3: " << moments1[4] << endl;
	cerr << "M00: " << moments1[0] << endl;
	cerr << "11.0/6.0: " << FracMom2(11.0/6.0, 0, moments1) << endl;
	cerr << "-1.0/2.0: " << FracMom2(-1.0/2.0, 0, moments1) << endl;
	cin.get();
}*/

	return sqrt(phi0 * phi1) * c;
}

Double TSoot::GetPhi( int x, int y, Double temp, Double *moments1, Double *moments2 )
{
	Double	phi0, phi1;
	Double	c = GetC( temp );
	
	phi0 = 		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 );

	phi1 = 		FracMom2( 7.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( 5.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( 5.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 7.0 / 6.0 + y, moments2 );

	return sqrt( phi0 * phi1 ) * c;
}

Double TSoot::GetPhiPAH( int x, int y, Double temp, Double *moments1, Double *moments2 )
{
	Double	phi0, phi1;
	Double	c = GetCPAH( temp );
	
	phi0 = 		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 );

	phi1 = 		FracMom2( 7.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( 5.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( 5.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 7.0 / 6.0 + y, moments2 );

	return sqrt( phi0 * phi1 ) * c;
}

//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i// 
Double TSoot::GetPhiPAHFROMA4( Double x, Double temp, Double *moments )
{
	Double	phi0, phi1;
	Double	c = GetC( temp );
	
	phi0 =      pow( 3.0, 1.0 / 3.0 ) * FracMom2(             x, moments )  
		+ 2.0 * pow( 3.0,-1.0 / 3.0 ) * FracMom2( 1.0 / 3.0 + x, moments ) 
		+       1.0 / 3.0             * FracMom2( 2.0 / 3.0 + x, moments );

	phi1 = 	    pow( 3.0, 7.0 / 3.0 ) * FracMom2(             x, moments )  
		+ 2.0 * pow( 3.0, 5.0 / 3.0 ) * FracMom2( 1.0 / 3.0 + x, moments ) 
		+       3.0                   * FracMom2( 2.0 / 3.0 + x, moments ) 
		+       pow( 3.0, 1.0 / 3.0 ) * FracMom2( 1.0       + x, moments )		
		+ 2.0 * pow( 3.0,-1.0 / 3.0 ) * FracMom2( 4.0 / 3.0 + x, moments )
		+       1.0 / 3.0             * FracMom2( 5.0 / 3.0 + x, moments );	

	return sqrt( phi0 * phi1 ) * c;
}
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//

Double TSoot::GetDerivPhiCond( int l, int x, int y, Double temp, Double *moments1, Double *moments2  )
{
	Double	phi0, phi1, phi, dPhi0dMl, dPhi1dMl;
	Double	c = GetC( temp );

	phi0 = 		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 );

	phi1 = 		FracMom2( 7.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( 5.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( 5.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 7.0 / 6.0 + y, moments2 );

	dPhi0dMl = 		GetDerivFracMom( l, 1.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * GetDerivFracMom( l, -1.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, -1.0 / 2.0 + x, moments1 ) 
					* FracMom2( 1.0 / 6.0 + y, moments2 );

	dPhi1dMl = 		GetDerivFracMom( l, 7.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * GetDerivFracMom( l, 5.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, 1.0 / 2.0 + x, moments1 ) 
					* FracMom2( 1.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, 1.0 / 6.0 + x, moments1 ) 
					* FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * GetDerivFracMom( l, -1.0 / 6.0 + x, moments1 ) 
					* FracMom2( 5.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, -1.0 / 2.0 + x, moments1 ) 
					* FracMom2( 7.0 / 6.0 + y, moments2 );

	phi = sqrt( phi0 * phi1 ) * c;
	
	if ( phi >= c * SMALLSOOT ) {
		return 0.5 * phi * ( dPhi0dMl / phi0 + dPhi1dMl / phi1 );
	}
	else {
		return 0.0;
	}
}

void TSoot::ComputePhiPAH( Double temp )
{
	Double	c = GetC( temp );
	Double	*fm = fFracMomPAH->vec;
	Double	*phi = fPhiPAH->vec;

	phi[k000] = ( fm[k01_6] * fm[km1_2] + fm[km1_6] * fm[km1_6] ) * 2.0;
	phi[k100] = ( fm[k07_6] * fm[km1_2] + fm[k05_6] * fm[km1_6] * 2.0 
									+ fm[k01_2] * fm[k01_6] ) * 2.0;
	
	phi[k011] = ( fm[k07_6] * fm[k01_2] + fm[k05_6] * fm[k05_6] ) * 2.0;
	phi[k111] = ( fm[k13_6] * fm[k01_2] + fm[k11_6] * fm[k05_6] * 2.0 
									+ fm[k03_2] * fm[k07_6] ) * 2.0;
	
	phi[k012] =   fm[k07_6] * fm[k03_2] + fm[k05_6] * fm[k11_6] * 2.0 
									+ fm[k01_2] * fm[k13_6];
	phi[k112] =   fm[k13_6] * fm[k03_2] 
					+ ( fm[k11_6] * fm[k11_6] + fm[k05_6] * fm[k17_6] ) * 2.0
				 	+ fm[k03_2] * fm[k13_6] + fm[k07_6] * fm[k05_2] + fm[k01_2] * fm[k19_6];
	
	phi[kh00] = sqrt( phi[k000] * phi[k100] ) * c;
	phi[kh11] = sqrt( phi[k011] * phi[k111] ) * c;
	phi[kh12] = sqrt( phi[k012] * phi[k112] ) * c;
	
/*	cerr << "phi[kh00] = " << phi[kh00] << TAB
		 << "phi[kh11] = " << phi[kh11] << TAB
		 << "phi[kh12] = " << phi[kh12] << NEWL;*/
}

Double TSoot::SourceCoagulation( int i )
{
	switch ( i ) {
		case 0:
			return -0.5 * fPhi->vec[kh00];
		case 1:
			return 0.0;
		case 2:
			return fPhi->vec[kh11];
		case 3:
			return 3.0 * fPhi->vec[kh12];
		default:
			cerr << "#error: no need to compute coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::SourceSurfDepCoag( int i, Double *mom, Double *Y, Double temp
					, Double density, Double *molarMass )
{
	Double	ksg, c, coag, S;
	Double	Cst = 1.0 / ( 6.0 * fCoagFact );

	switch ( i ) {
		case 0:
			ksg = GetSurfGrowthCoeff( Y, density, molarMass );
			c = GetC( temp );
			coag = -0.5 * fPhi->vec[kh00];
			S = MAX( 1.0e-30, GetSCoag( mom, temp ) );
//			fprintf( stderr, "S = %g\tCoag = %g\tksg = %g\tc = %g\tCst = %g\n", S, coag, ksg, c, Cst );

			return coag * ( 1.0 - 1.0 / ( 1.0 + Cst * fPhi->vec[kh00] * ksg / ( c * c * S ) ) );
		case 1:
			return 0.0;
		case 2:
			fprintf( stderr, "###error: SourceSurfDepCoag for higher moments not yet implemented\n" );
			return 0.0;
		case 3:
			fprintf( stderr, "###error: SourceSurfDepCoag for higher moments not yet implemented\n" );
			return 0.0;
		default:
			cerr << "#error: no need to compute coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::GetSCoag( Double *mom, Double /*temp*/ )
{
	Double	S00, S01, S10, S11;
	Double	M_1_6 = FracMom2( 1.0 / 6.0, mom );
	Double	M_m1_6 = FracMom2( -1.0 / 6.0, mom );
	Double	M_1_2 = FracMom2( 1.0 / 2.0, mom );
	Double	M_m1_2 = FracMom2( -1.0 / 2.0, mom );
	Double	M_3_2 = FracMom2( 3.0 / 2.0, mom );
	Double	M_5_2 = FracMom2( 5.0 / 2.0, mom );
	Double	M_5_6 = FracMom2( 5.0 / 6.0, mom );
	Double	M_11_6 = FracMom2( 11.0 / 6.0, mom );
	Double	M_17_6 = FracMom2( 17.0 / 6.0, mom );
	
	S00 = M_m1_6 * M_m1_6 * ( M_m1_2 + 2.0 * M_1_2 + M_3_2 );
	S01 = 2.0 * ( M_11_6 * M_m1_6 * M_m1_2 
					+ M_5_6 * M_5_6 * M_m1_2 
					+ 2.0 * M_5_6 * M_m1_6 * M_1_2 ) 
			+ M_m1_6 * M_m1_6 * M_3_2;
	S10 = 2.0 * M_5_6 * M_m1_6 * M_m1_2
			+ M_m1_6 * M_m1_6 * M_1_2
			+ 4.0 * M_5_6 * M_m1_6 * M_1_2
			+ 2.0 * M_m1_6 * M_m1_6 * M_3_2
			+ 2.0 * M_5_6 * M_m1_6 * M_3_2
			+ M_m1_6 * M_m1_6 * M_5_2;
	S11 = 2.0 * M_17_6 * M_m1_6 * M_m1_2
			+ M_m1_6 * M_m1_6 * M_5_2
			+ 6.0 * M_11_6 * M_5_6 * M_m1_2
			+ 6.0 * M_11_6 * M_m1_6 * M_1_2
			+ 6.0 * M_5_6 * M_m1_6 * M_3_2
			+ 6.0 * M_5_6 * M_5_6 * M_1_2;


	return pow( S00 * S10, 1.0/3.0 ) * pow( S01 * S11, 1.0/6.0 );
}

//Mueller (9/7/07)
#ifdef VS
Double TSoot::SourceCoagulationNew(int k, Double temp, Double *moments)
{
	Double	FitExp = -0.1848;
	Double	FitV = -2.0 * FitExp;
	Double	FitS = 3.0 * FitExp;
	Double	FitC = 0.6324;

/*	double complex = 0.5 * GetPhi(temp, moments, moments, 3.0/2.0) - GetPhi(1, 0, 0, 0, temp, 3.0, moments, moments);
	double simple = 0.0;
	if (simple < 0.999*complex || simple > 1.001*complex)
	{
		cerr << "Complex term: " << complex << endl;
		cerr << "Simple term: " << simple << endl;
		cin.get();
	}*/

	switch (k)
	{
		case 0:	//M(0,0)
//			return -0.5 * GetPhi(0, 0, 0, 0, temp, 1.8, moments, moments); //Aggregation, d(dc)=d(dp)
			return -0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments); //Coalescence
		case 1:	//M(1,0)
			return 0.0; //Aggregation, d(dc)=d(dp), Coalescence
		case 2:	//M(0,1)
//			return 0.0; //Aggregation
//			return 0.5 * FitC * GetPhi(FitV-1.0, FitS+1.0, 1, 0, temp, 1.8, moments, moments) - 0.5 * GetPhi(0, 1, 0, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
			return 0.5 * GetPhi(temp, moments, moments, 7.0/6.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence
		case 3: //M(-1/2,0)
//			return 0.5 * GetPhi(temp, moments, moments, 0) - GetPhi(-1.0/2.0, 0, 0, 0, temp, 3.0, moments, moments); //Coalescence
//			return GetPhi(1, 0, 1, 0, temp, 3.0, moments, moments);
			return 0.5 * GetPhi(temp, moments, moments, 5.0/6.0) - GetPhi(0, 1.0/2.0, 0, 0, temp, 3.0, moments, moments);
		case 4: //M(0,2)
//			return 0.5 * GetPhi(temp, moments, moments, 11.0/6.0) - GetPhi(0, 2, 0, 0, temp, 3.0, moments, moments); //Coalescence
			return 0.5 * GetPhi(temp, moments, moments, 5.0/6.0) - GetPhi(0, 1.0/2.0, 0, 0, temp, 3.0, moments, moments);
		case 5: //M(7/6,1)
//			return 0.5 * GetPhi(temp, moments, moments, 7.0/3.0) - GetPhi(7.0/6.0, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence
//			return 3.0 * GetPhi(2, 0, 1, 0, temp, 3.0, moments, moments);
			return 0.5 * GetPhi(temp, moments, moments, 7.0/6.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments);
		default:
			cerr << "#error: coagulation not yet implemented for " << k << " moments" << endl;
			exit(2);
	}

	return 0.0;
}
#endif

//Mueller (7/20/07)
Double TSoot::SourceCoagulationNew(int i, int j, Double temp, Double *moments)
{
	Double	FitExp = -0.1848;
	Double	FitV = -2.0 * FitExp;
	Double	FitS = 3.0 * FitExp;
	Double	FitC = 0.6324;

/*	if (moments[1] / moments[0] < 1)
	moments[0] /= 1000000;

	double first = 0.5 * GetPhi(temp, moments, moments, (double(i)/6.0)+(2.0/3.0)*(double(j)/6.0)+(1.0/2.0));
	double second = GetPhi(double(i)/6.0, double(j)/6.0, 0, 0, temp, 3.0, moments, moments);
	double diff = first - second;

	cerr << "First Term: " << first << "\tSecond Term: " << second << "\tDifference: " << diff << endl;
	if (isnan(first) || isnan(second) || (diff > 0 && i==0 && j==0))
	{
		cerr << "Moment(0,0): " << moments[0] << endl;
		cerr << "Moment(1,0): " << moments[1] << endl;
		cerr << "Moment(0,1): " << moments[2] << endl;
		cerr << "Moment(2,0): " << moments[3] << endl;
		cerr << "Moment(0,2): " << moments[4] << endl;
		cerr << "Moment(3,0): " << moments[5] << endl;
		cin.get();
	}*/

	return 2.0 * (0.5 * GetPhi(temp, moments, moments, (double(i)/6.0)+(2.0/3.0)*(double(j)/6.0)+(1.0/2.0)) - GetPhi(double(i)/6.0, double(j)/6.0, 0, 0, temp, 3.0, moments, moments));

	if (i==0 && j==0)
	{
//		return -0.5 * GetPhi(0, 0, 0, 0, temp, 1.8, moments, moments); //Aggregation, d(dc)=d(dp)
		return -0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==6 && j==0)
	{
		return 0.0;
	}
	else if (i==0 && j==6)
	{
//		return 0.0; //Aggregation
//		return 0.5 * FitC * GetPhi(FitV-1.0, FitS+1.0, 1, 0, temp, 1.8, moments, moments) - 0.5 * GetPhi(0, 1, 0, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
		return 0.5 * GetPhi(temp, moments, moments, 7.0/6.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==0 && j==3)
		return 0.5 * GetPhi(temp, moments, moments, 5.0/6.0) - GetPhi(0, 1.0/2.0, 0, 0, temp, 3.0, moments, moments);
	else if (i==2 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, 5.0/6.0) - GetPhi(1.0/3.0, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==4 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, 7.0/6.0) - GetPhi(2.0/3.0, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==2 && j==3)
		return 0.5 * GetPhi(temp, moments, moments, 7.0/6.0) - GetPhi(1.0/3.0, 1.0/2.0, 0, 0, temp, 3.0, moments, moments);
	else if (i==6 && j==6)
		return 0.5 * GetPhi(temp, moments, moments, 13.0/6.0) - GetPhi(1, 1, 0, 0, temp, 3.0, moments, moments);
	else if (i==12 && j==0)
		return GetPhi(1, 0, 1, 0, temp, 3.0, moments, moments);
/*	else if (i==-12 && j==18)
		return -0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments);// - GetPhi(-2, 3, 0, 0, temp, 3.0, moments, moments);*/
	else if (i==11 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, 7.0/3.0) - GetPhi(11.0/6.0, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==7 && j==6)
		return 0.5 * GetPhi(temp, moments, moments, 7.0/3.0) - GetPhi(7.0/6.0, 1, 0, 0, temp, 3.0, moments, moments);
	else if (i==0 && j==18)
		return GetPhi(1, 0, 1, 0, temp, 3.0, moments, moments);
	else if (i==-3 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, 0) - GetPhi(-1.0/2.0, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==13 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, 8.0/3.0) - GetPhi(13.0/6.0, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==7 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, 5.0/3.0) - GetPhi(7.0/6.0, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==4 && j==3)
		return 0.5 * GetPhi(temp, moments, moments, 3.0/2.0) - GetPhi(2.0/3.0, 1.0/2.0, 0, 0, temp, 3.0, moments, moments);
//		return GetPhi(1, 0, 0, 0, temp, 3.0, moments, moments) - GetPhi(2.0/3.0, 1.0/2.0, 0, 0, temp, 3.0, moments, moments);
	else if (i==0 && j==-6)
		return 0.5 * GetPhi(temp, moments, moments, -1.0/6.0) - GetPhi(0, -1, 0, 0, temp, 3.0, moments, moments);
	else if (i==6 && j==-6)
		return 0.5 * GetPhi(temp, moments, moments, 5.0/6.0) - GetPhi(1, -1, 0, 0, temp, 3.0, moments, moments);
	else if (i==12 && j==-6)
		return 0.5 * GetPhi(temp, moments, moments, 11.0/6.0) - GetPhi(2, -1, 0, 0, temp, 3.0, moments, moments);
	else if (i==-6 && j==-6)
		return 0.5 * GetPhi(temp, moments, moments, -7.0/6.0) - GetPhi(-1, -1, 0, 0, temp, 3.0, moments, moments);
	else if (i==-6 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, -1.0/2.0) - GetPhi(-1, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==-6 && j==6)
		return 0.5 * GetPhi(temp, moments, moments, 1.0/6.0) - GetPhi(-1, 1, 0, 0, temp, 3.0, moments, moments);
	else if (i==-6 && j==12)
		return 0.5 * GetPhi(temp, moments, moments, 5.0/6.0) - GetPhi(-1, 2, 0, 0, temp, 3.0, moments, moments);
	else if (i==-12 && j==-6)
		return 0.5 * GetPhi(temp, moments, moments, -13.0/6.0) - GetPhi(-2, -1, 0, 0, temp, 3.0, moments, moments);
	else if (i==-12 && j==0)
		return 0.5 * GetPhi(temp, moments, moments, -3.0/2.0) - GetPhi(-2, 0, 0, 0, temp, 3.0, moments, moments);
	else if (i==-12 && j==6)
		return 0.5 * GetPhi(temp, moments, moments, -5.0/6.0) - GetPhi(-2, 1, 0, 0, temp, 3.0, moments, moments);
	else if (i==-12 && j==12)
		return 0.5 * GetPhi(temp, moments, moments, -1.0/6.0) - GetPhi(-2, 2, 0, 0, temp, 3.0, moments, moments);
	else if (i==-12 && j==18)
	{
/*		cerr << "First Term: " << 0.5 * GetPhi(temp, moments, moments, 1.0/2.0) << endl;
		cerr << "Second Term: " << GetPhi(-2, 3, 0, 0, temp, 3.0, moments, moments);
		if (isnan(GetPhi(-2, 3, 0, 0, temp, 3.0, moments, moments)) || isinf(GetPhi(-2, 3, 0, 0, temp, 3.0, moments, moments)))
		cin.get();*/
		return 0.5 * GetPhi(temp, moments, moments, 1.0/2.0) - GetPhi(-2, 3, 0, 0, temp, 3.0, moments, moments);
	}
	else if (i==0 && j==12)
		return 0.5 * GetPhi(temp, moments, moments, 11.0/6.0) - GetPhi(0, 2, 0, 0, temp, 3.0, moments, moments);
/*	else if (i==2 && j==0)
	{
//		return GetPhi(1, 0, 1, 0, temp, 1.8, moments, moments); //Aggregation, d(dc)=d(dp)
		return GetPhi(1, 0, 1, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==0 && j==2)
	{
//		return GetPhi(0, 1, 0, 1, temp, 1.8, moments, moments); //Aggregation
//		return FitC * GetPhi(FitV-1.0, FitS+2.0, 1, 0, temp, 1.8, moments, moments) + 0.5 * FitC * FitC * GetPhi(2.0*FitV-2.0, 2.0*FitS+2.0, 2, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
		return 0.5 * GetPhi(temp, moments, moments, 11.0/6.0) - GetPhi(0, 2, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==1 && j==1)
	{
//		return GetPhi(1, 0, 0, 1, temp, 1.8, moments, moments); //Aggregation
//		return 0.5 * FitC * GetPhi(FitV, FitS+1.0, 1, 0, temp, 1.8, moments, moments) + 0.5 * GetPhi(0, 1, 1, 0, temp, 1.8, moments, moments) + 0.5 * FitC * GetPhi(FitV-1.0, FistS+1.0, 2, 0, temp, 1.8, moments, moments) - 0.5 * GetPhi(1, 1, 0, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
		return 0.5 * GetPhi(temp, moments, moments, 13.0/6.0) - GetPhi(1, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==3 && j==0)
	{
//		return 3.0 * GetPhi(2, 0, 1, 0, temp, 1.8, moments, moments); //Aggregation, d(dc)=d(dp)
		return 3.0 * GetPhi(2, 0, 1, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==0 && j==3)
	{
//		return 3.0 * GetPhi(0, 2, 0, 1, temp, 1.8, moments, moments); //Aggregation
//		return 1.5 * FitConst * GetPhi(FitV-1.0, FitS+3.0, 1, 0, temp, 1.8, moments, moments) + 1.5 * FitC * FitC * GetPhi(2.0*FitV-2.0, 2.0*FitS+3.0, 2, 0, temp, 1.8, moments, moments) + 0.5 * FitC * FitC * FitC * GetPhi(3.0*FitV-3.0, 3.0*FitS+3.0, 3, 0, temp, 1.8, moments, moments) - 0.5 * GetPhi(0, 3, 0, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
		return GetPhi(2, 0, 0, 0, temp, 3.0, moments, moments) + GetPhi(1, 0, 1, 0, temp, 3.0, moments, moments) - GetPhi(0, 3, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==2 && j==1)
	{
//		return GetPhi(2, 0, 0, 1, temp, 1.8, moments, moments) + 2.0 * GetPhi(1, 1, 1, 0, temp, 1.8, moments, moments); //Aggregation
//		return 0.5 * FitC * GetPhi(FitV+1.0, FitS+1.0, 1, 0, temp, 1.8, moments, moments) + GetPhi(1, 1, 1, 0, temp, 1.8, moments, moments) + FitC * GetPhi(FitV, FitS+1.0, 2, 0, temp, 1.8, moments, moments) + 0.5 * GetPhi(0, 1, 2, 0, temp, 1.8, moments, moments) + 0.5 * FitC * GetPhi(FitV-1.0, FitS+1.0, 3, 0, temp, 1.8, moments, moments) - 0.5 * GetPhi(2, 1, 0, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
		return 0.5 * GetPhi(temp, moments, moments, 19.0/6.0) - GetPhi(2, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==1 && j==2)
	{
//		return GetPhi(1, 0, 0, 2, temp, 1.8, moments, moments) + 2.0 * GetPhi(1, 1, 0, 1, temp, 1.8, moments, moments); //Aggregation
//		return FitC * GetPhi(FitV, FitS+2.0, 1, 0, temp, 1.8, moments, moments) + 0.5 * FitC * FitC * GetPhi(2.0*FitV-1.0, 2.0*FitS+2.0, 2, 0, temp, 1.8, moments, moments) + 0.5 * GetPhi(0, 2, 1, 0, temp, 1.8, moments, moments) + FitC * GetPhi(FitV-1.0, FitS+2.0, 2, 0, temp, 1.8, moments, moments) + 0.5 * FitC * FitC * GetPhi(2.0*FitV-2.0, 2.0*FitS+2.0, 3, 0, temp, 1.8, moments, moments) - 0.5 * GetPhi(1, 2, 0, 0, temp, 1.8, moments, moments); //d(dc)=d(dp)
		return 0.5 * GetPhi(temp, moments, moments, 17.0/6.0) - GetPhi(1, 2, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}*/
	else
	{
		cerr << "#error: no need to compute source term for coagulation for M_" << i << "," << j << endl;
		exit(2);
	}

	return 0.0;
}

#ifndef VS
Double TSoot::SourceCoagulationNew( int i, Double temp, Double *moments )
{
	switch ( i ) {
		case 0:
//			return ( moments[0] > 1.0e-30 ) ? -0.5 * GetPhi( 0, 0, temp, moments ) : 0.0;
			return -/* 1.0 / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[0], 1.0e-60 ) ) * */ 0.5 * GetPhi( 0, 0, temp, moments );
		case 1:
			return 0.0;
		case 2:
			return GetPhi( 1, 1, temp, moments );
		case 3:
			return 3.0 * GetPhi( 1, 2, temp, moments );
		default:
			cerr << "#error: no need to compute coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}
#endif

Double TSoot::SourceCondensation( int i, Double temp, Double *pahMoments, Double *moments )
{
	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return GetBetaCond( temp, moments ) * moments[0] * pahMoments[1];
		case 2:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[2] 
								+ 2.0 * moments[1] * pahMoments[1] );
		case 3:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[3] 
								+ 3.0 * ( moments[1] * pahMoments[2] 
										+ moments[2] * pahMoments[1] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

//Mueller (9/7/07)
/*Double TSoot::SourceCondensationNew(int k, Double temp, Double *pahMoments, Double *moments)
{
	Double	FitExp = -0.1848;
	Double	FitV = -2.0 * FitExp;
	Double	FitS = 3.0 * FitExp;
	Double	FitC = 0.6324;

	switch (k)
	{
		case 0:	//M(0,0)
			return 0.0; //d(dc)=d(dp), d(np)=0, Coalescence
		case 1:	//M(1,0)
//			return GetPhi(0, 0, 1, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp), d(np)=0
			return GetPhi(0, 0, 1, 0, temp, 3.0, moments, pahMoments); //Coalescence
		case 2:	//M(0,1)
//			retur2.0 * GetPhi(1, 0, 1, 0, temp, 3.0, moments, pahMoments)n FitC * GetPhi(FitV-1.0, FitS+1.0, 1, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//			return (2.0/3.0) * GetPhi(-1, 1, 1, 0, temp, 1.8, moments, pahMoments); //d(np)=0
			return GetPhi(temp, moments, pahMoments, 7.0/6.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
		case 3:	//M(-1/2,0)
			return GetPhi(temp, moments, pahMoments, 0) - GetPhi(-1.0/2.0, 0, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
		case 4:	//M(0,2)
			return GetPhi(temp, moments, pahMoments, 11.0/6.0) - GetPhi(0, 2, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
		case 5:	//M(7/6,1)
			return GetPhi(temp, moments, pahMoments, 7.0/3.0) - GetPhi(7.0/6.0, 1, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
		default:
			cerr << "#error: condensation not yet implemented for " << k << " moments" << endl;
			exit(2);

	return 0.0;
}*/

//Mueller--Updated 9/7/07
Double TSoot::SourceCondensationNew(int i, int j, Double temp, Double *pahMoments, Double *moments)
{	
	Double	FitExp = -0.1848;
	Double	FitV = -2.0 * FitExp;
	Double	FitS = 3.0 * FitExp;
	Double	FitC = 0.6324;

	Double ii = double(i)/6.0;
	Double jj = double(j)/6.0;

/*	double complex = 2.0 * GetPhi(1, 0, 1, 0, temp, 1.8, moments, pahMoments) + GetPhi(0, 0, 2, 0, temp, 1.8, moments, pahMoments);
	double simple = 2.0 * GetPhi(1, 0, 1, 0, temp, 1.8, moments, pahMoments);
	if (simple < 0.9*complex || simple > 1.1*complex)
	{
		cerr << "Complex term: " << complex << endl;
		cerr << "Simple term: " << simple << endl;
		cin.get();
	}*/
	//Some difference with the taylor expansion, but is that necessarily worse?

	//Pure Coalescence
	return (ii+(2.0/3.0)*jj) * GetPhi(ii-1.0, jj, 1, 0, temp, 3.0, moments, pahMoments);

/*	cerr << "First term: " << 2.0 * GetPhi(1, 0, 1, 0, temp, 3.0, moments, pahMoments) << endl;
	cerr << "Second term: " << GetPhi(0, 0, 2, 0, temp, 3.0, moments, pahMoments) << endl;
	cin.get();*/
	
	if (i==0 && j==0)
	{
		return 0.0;
	}
	else if (i==6 && j==0)
	{
//		return GetPhi(0, 0, 1, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp), d(np)=0
		return GetPhi(0, 0, 1, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==0 && j==6)
	{
//		return FitC * GetPhi(FitV-1.0, FitS+1.0, 1, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//		return (2.0/3.0) * GetPhi(-1, 1, 1, 0, temp, 1.8, moments, pahMoments); //d(np)=0
		return GetPhi(temp, moments, pahMoments, 7.0/6.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==2 && j==3)
		return GetPhi(temp, moments, pahMoments, 7.0/6.0) - GetPhi(1.0/3.0, 1.0/2.0, 0, 0, temp, 3.0, moments, pahMoments);
	else if (i==7 && j==6)
		return GetPhi(temp, moments, pahMoments, 7.0/3.0) - GetPhi(7.0/6.0, 1, 0, 0, temp, 3.0, moments, pahMoments);
	else if (i==-3 && j==0)
		return GetPhi(temp, moments, pahMoments, 0) - GetPhi(-1.0/2.0, 0, 0, 0, temp, 3.0, moments, pahMoments);
	else if (i==13 && j==0)
		return GetPhi(temp, moments, pahMoments, 8.0/3.0) - GetPhi(13.0/6.0, 0, 0, 0, temp, 3.0, moments, pahMoments);
	else if (i==7 && j==0)
		return GetPhi(temp, moments, pahMoments, 5.0/3.0) - GetPhi(7.0/6.0, 0, 0, 0, temp, 3.0, moments, pahMoments);
	else if (i==4 && j==3)
		return GetPhi(temp, moments, pahMoments, 3.0/2.0) - GetPhi(2.0/3.0, 1.0/2.0, 0, 0, temp, 3.0, moments, pahMoments);
//		return GetPhi(1, 0, 0, 0, temp, 3.0, moments, pahMoments) + GetPhi(0, 0, 1, 0, temp, 3.0, moments, pahMoments) - GetPhi(2.0/3.0, 1.0/2.0, 0, 0, temp, 3.0, moments, pahMoments);
	else if (i==0 && j==12)
		return GetPhi(temp, moments, pahMoments, 11.0/6.0) - GetPhi(0, 2, 0, 0, temp, 3.0, moments, pahMoments);
/*	else if (i==2 && j==0)
	{
//		return 2.0 * GetPhi(1, 0, 1, 0, temp, 1.8, moments, pahMoments) + GetPhi(0, 0, 2, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp), d(np)=0
		return 2.0 * GetPhi(1, 0, 1, 0, temp, 3.0, moments, pahMoments) + GetPhi(0, 0, 2, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==0 && j==2)
	{
//		return 2.0 * FitC * GetPhi(FitV-1.0, FitS+2.0, 1, 0, temp, 1.8, moments, pahMoments) + FitC * FitC * GetPhi(2.0*FitV-2.0, 2.0*FitS+2.0, 2, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//		return (4.0/3.0) * GetPhi(-1, 2, 1, 0, temp, 1.8, moments, pahMoments) + (4.0/9.0) * GetPhi(-2, 2, 2, 0, temp, 1.8, moments, pahMoments); //d(np)=0
		return GetPhi(temp, moments, pahMoments, 11.0/6.0) - GetPhi(0, 2, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==1 && j==1)
	{
//		return FitC * GetPhi(FitV, FitS+1.0, 1, 0, temp, 1.8, moments, pahMoments) + GetPhi(0, 1, 1, 0, temp, 1.8, moments, pahMoments) + FitC * GetPhi(FitV-1.0, FitS+1.0, 2, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//		return (5.0/3.0) * GetPhi(0, 1, 1, 0, temp, 1.8, moments, pahMoments) + (2.0/3.0) * GetPhi(-1, 1, 2, 0, temp, 1.8, moments, pahMoments); //d(np)=0
		return GetPhi(temp, moments, pahMoments, 13.0/6.0) - GetPhi(1, 1, 0, 0, temp, 3.0, moments, pahMoments); //Colescence
	}
	else if (i==3 && j==0)
	{
//		return 3.0 * GetPhi(2, 0, 1, 0, temp, 1.8, moments, pahMoments) + 3.0 * GetPhi(1, 0, 2. 0, temp, 1.8, moments, pahMoments) + GetPhi(0, 0, 3, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp), d(np)=0
		return 3.0 * GetPhi(2, 0, 1, 0, temp, 3.0, moments, pahMoments) + 3.0 * GetPhi(1, 0, 2, 0, temp, 3.0, moments, pahMoments) + GetPhi(0, 0, 3, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==0 && j==3)
	{
//		return 3.0 * FitC * GetPhi(FitV-1.0, FitS+3.0, 1, 0, temp, 1.8, moments, pahMoments) + 3.0 * FitC * FitC * GetPhi(2.0*FitV-2.0, 2.0*FitS+3.0, 2, 0, temp, 1.8, moments, pahMoments) + FitC * FitC * FitC * GetPhi(3.0*FitV-3.0, 3.0*FitS+3.0, 3, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//		return 2.0 * GetPhi(-1, 3, 1, 0, temp, 1.8, moments, pahMoments) + (4.0/3.0) * GetPhi(-2, 3, 2, 0, temp, 1.8, moments, pahMoments) + (8.0/27.0) * GetPhi(-3, 3, 3, 0, temp, 1.8, moments, pahMoments); //d(np)=0
		return GetPhi(2, 0, 0, 0, temp, 3.0, moments, pahMoments) + 2.0 * GetPhi(1, 0, 1, 0, temp, 3.0, moments, pahMoments) + GetPhi(0, 0, 2, 0, temp, 3.0, moments, pahMoments) - GetPhi(0, 3, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==2 && j==1)
	{
//		return  FitC * GetPhi(FitV+1.0, FitS+1.0, 1, 0, temp, 1.8, moments, pahMoments) + 2.0 * GetPhi(1, 1, 1, 0, temp, 1.8, moments, pahMoments) + 2.0 * FitC * GetPhi(FitV, FitS+1.0, 2, 0, temp, 1.8, moments, pahMoments) + GetPhi(0, 1, 2, 0, temp, 1.8, moments, pahMoments) + FitC * GetPhi(FitV-1.0, FitS+1.0, 3, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//		return (8.0/3.0) * GetPhi(1, 1, 1, 0, temp, 1.8, moments, pahMoments) + (7.0/3.0) * GetPhi(0, 1, 2, 0, temp, 1.8, moments, pahMoments) + (2.0/3.0) * GetPhi(-1, 1, 3, 0, temp, 1.8, moments, pahMoments); //d(np)=0
		return GetPhi(temp, moments, pahMoments, 19.0/6.0) - GetPhi(2, 1, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}
	else if (i==1 && j==2)
	{
//		return 2.0 * FitC * GetPhi(FitV, FitS+2.0, 1, 0, temp, 1.8, moments, pahMoments) + FitC * FitC * GetPhi(2.0*FitV-1.0, 2.0*FitS+2.0, 2, 0, temp, 1.8, moments, pahMoments) + GetPhi(0, 2, 1, 0, temp, 1.8, moments, pahMoments) + 2.0 * FitC * GetPhi(FitV-1.0, FitS+2.0, 2, 0, temp, 1.8, moments, pahMoments) + FitC * FitC * GetPhi(2.0*FitV-2.0, 2.0*FitS+2.0, 3, 0, temp, 1.8, moments, pahMoments); //d(dc)=d(dp)
//		return (7.0/3.0) * GetPhi(0, 2, 1, 0, temp, 1.8, moments, pahMoments) + (16.0/9.0) * GetPhi(-1, 2, 2, 0, temp, 1.8, moments, pahMoments) + (4.0/9.0) * GetPhi(-2, 2, 3, 0, temp, 1.8, moments, pahMoments); //d(np)=0
		return GetPhi(temp, moments, pahMoments, 17.0/6.0) - GetPhi(1, 2, 0, 0, temp, 3.0, moments, pahMoments); //Coalescence
	}*/
	else
	{
		cerr << "#error: no need to compute source term for condensation for M_" << i << "," << j << endl;
		exit( 2 );
	}

	return 0.0;
}

#ifdef NEWPOLYNUCCOND
Double TSoot::SourceCondensationNew( int i, Double temp, Double *pahMoments, Double *moments, Double */*Y*/
                                   , Double /*density*/, Double */*molarMass*/  )
{

	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return GetBetaCond( temp, moments ) * moments[0] * pahMoments[1];
		case 2:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[2] 
								+ 2.0 * moments[1] * pahMoments[1] );
		case 3:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[3] 
								+ 3.0 * ( moments[1] * pahMoments[2] 
										+ moments[2] * pahMoments[1] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
#else
Double TSoot::SourceCondensationNew( int i, Double temp, Double *pahMom, Double *mom, Double *Y
                                   , Double density, Double *molarMass  )
{	
	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
#	ifdef PAHFROMA4
// changed by hp from
//			return 9.0 * pahMom[0] * GetPhiPAHFROMA4( x, temp, pahMom );
// to
			return 9.0 * pahMom[0] * GetPhiPAHFROMA4( -0.5, temp, mom );
#	else
			return GetPhi( 0, 1, temp, mom, pahMom );
#	endif
		case 2:
			return GetPhi( 0, 2, temp, mom, pahMom ) + 2.0 * GetPhi( 1, 1, temp, mom, pahMom );
		case 3:
			return GetPhi( 0, 3, temp, mom, pahMom ) 
				+ 3.0 * ( GetPhi( 1, 2, temp, mom, pahMom ) + GetPhi( 2, 1, temp, mom, pahMom ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
#endif
}

Double TSoot::SourceSurfGrowth( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSurfGrowthCoeff( Y, density, molarMass );
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return coeff * fracMom[k23_24];
		case 2:
			return coeff * ( fracMom[k23_24] + 2.0 * fracMom[k47_24] );
		case 3:
			return coeff * ( fracMom[k23_24] + 3.0 * ( fracMom[k47_24] + fracMom[k71_24] ) );
		default:
		cerr << "#error: no need to compute source term for surface growth for M_" << i << NEWL;
			exit( 2 );

	}
	return 0.0;
}

//Mueller (8/20/07)
Double TSoot::GetThetaSG(Double frac, Double *moments)
{
	Double	theta0, theta1, theta2, theta3;

	theta0 = moments[2];
	theta1 = FracMom2(1.0, 1.0, moments) + moments[2];

	if (fNSootMoments > 3)
		theta2 = FracMom2(2.0, 1.0, moments) + 2.0 * FracMom2(1.0, 1.0, moments) + moments[2];

	if (fNSootMoments > 6)
		theta3 = FracMom2(3.0, 1.0, moments) + 3.0 * FracMom2(2.0, 1.0, moments) + 3.0 * FracMom2(1.0, 1.0, moments) + moments[2];

	if (fNSootMoments==3 || fNSootMoments==6)
		return pow(theta0, 1.0-frac) * pow(theta1, frac);
	else if (fNSootMoments==6)
		return pow(theta0, (1.0-frac)*(1.0-0.5*frac)) * pow(theta1, frac*(2.0-frac)) * pow(theta2, 0.5*frac*(frac-1.0));
	else if (fNSootMoments==10)
		return pow(theta0, (1.0-frac)*(1.0-0.5*frac)*(1.0-(1.0/3.0)*frac)) * pow(theta1, frac*(2.0-frac)*(1.5-0.5*frac)) * pow(theta2, 0.5*frac*(frac-1.0)*(3.0-frac)) * pow(theta3, (1.0/3.0)*frac*(0.5*frac-0.5)*(frac-2.0));

	return 0.0;
}

//Mueller (7/20/07)
Double TSoot::SourceSurfGrowthNew(int i, int j, Double *moments, Double *Y, Double density, Double *molarMass)
{
	Double	coeffF = GetSurfGrowthCoeffForw(Y, density, molarMass);
	Double	coeffB = GetSurfGrowthCoeffBackw(Y, density, molarMass);

	Double	ChiPrime = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

	Double	FitExp = -0.1848;
	Double	FitV = -2.0 * FitExp;
	Double	FitS = 3.0 * FitExp;
	Double	FitC = 0.6324;

	double complex = (coeffF-coeffB) * ChiPrime * (2.0 * FracMom2(1.0, 1.0, moments) + moments[2]);
	double simple = (coeffF-coeffB) * ChiPrime * (2.0 * FracMom2(1.0, 1.0, moments));
	if (simple < 0.999*complex || simple > 1.001*complex)
	{
		cerr << "Complex term: " << complex << endl;
		cerr << "Simple term: " << simple << endl;
		cin.get();
	}

	if (i==0 && j==0)
	{
//		return 0.0;
		return -coeffB * ChiPrime * FracMom2(-1.0, 1.0, moments);
	}
	else if (i==6 && j==0)
	{
		return (coeffF-coeffB) * ChiPrime * moments[2];
	}
	else if (i==0 && j==6)
	{
//		return FitC * (coeffF-coeffB) * ChiPrime * FracMom2(FitV-1.0, FitS+2.0, moments); //d(dc)=d(dp)
//		return (2.0/3.0) * (coeffF-coeffB) * ChiPrime * FracMom2(-1.0, 2.0, moments); //d(np)=0
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(2.0/3.0, moments) - FracMom2(0.0, 2.0, moments)); //Coalescence
	}
	else if (i==2 && j==3)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(2.0/3.0, moments) - FracMom2(1.0/3.0, 3.0/2.0, moments));
	else if (i==7 && j==6)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(11.0/6.0, moments) - FracMom2(7.0/6.0, 2.0, moments));
	else if (i==-3 && j==0)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(-1.0/2.0, moments) - FracMom2(-1.0/2.0, 1.0, moments));
	else if (i==13 && j==0)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(13.0/6.0, moments) - FracMom2(13.0/6.0, 1.0, moments));
	else if (i==7 && j==0)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(7.0/6.0, moments) - FracMom2(7.0/6.0, 1.0, moments));
	else if (i==4 && j==3)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(1.0, moments) - FracMom2(2.0/3.0, 3.0/2.0, moments));
	else if (i==0 && j==12)
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(4.0/3.0, moments) - FracMom2(0.0, 3.0, moments));
/*	else if (i==2 && j==0)
	{
		return (coeffF-coeffB) * ChiPrime * (2.0 * moments[4] + moments[1]);
	}
	else if (i==0 && j==2)
	{
//		return (coeffF-coeffB) * ChiPrime * (2.0 * FitC * FracMom2(FitV-1.0, FitS+3.0, moments) + FitC * FitC * FracMom2(2.0*FitV-2.0, 2.0*FitS+3.0, moments)); //d(dc)=d(dp)
//		return (coeffF-coeffB) * ChiPrime * ((4.0/3.0) * FracMom2(-1.0, 3.0, moments) + (4.0/9.0) * FracMom2(-2.0, 3.0, moments)); //d(np)=0
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(4.0/3.0, moments) - FracMom2(0.0, 3.0, moments)); //Coalescence
	}
	else if (i==1 && j==1)
	{
//		return (coeffF-coeffB) * ChiPrime * (FitC * FracMom2(FitV, FitS+2.0, moments) + moments[3] + FitC * FracMom2(FitV-1.0, FitS+2.0, moments)); //d(dc)=d(dp)
//		return (coeffF-coeffB) * ChiPrime * ((5.0/3.0) * moments[3] + (2.0/3.0) * FracMom2(-1.0, 2.0, moments)); //d(np)=0
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(5.0/3.0, moments) - FracMom2(1.0, 2.0, moments)); //Coalescence
	}
	else if (i==3 && j==0)
	{
		return (coeffF-coeffB) * ChiPrime * (3.0 * moments[8] + 3.0 * moments[4] + moments[1]);
	}
	else if (i==0 && j==3)
	{
//		return (coeffF-coeffB) * ChiPrime * (3.0 * FitC * FracMom2(FitV-1.0, FitS+4.0, moments) + 3.0 * FitC * FitC * FracMom2(2.0*FitV-2.0, 2.0*FitS+4.0, moments) + FitC * FitC * FitC * FracMom2(3.0*FitV-3.0, 3.0*FitS+4.0, moments)); //d(dc)=d(dp)
//		return (coeffF-coeffB) * ChiPrime * (2.0 * FracMom2(-1.0, 4.0, moments) + (4.0/3.0) * FracMom2(-2.0, 4.0, moments) + (8.0/27.0) * FracMom2(-3.0, 4.0, moments)); //d(np)=0
		return (coeffF-coeffB) * ChiPrime * (moments[8] + 2.0 * moments[4] + moments[1] - FracMom2(0.0, 4.0, moments)); //Coalescence
	}
	else if (i==2 && j==1)
	{
//		return (coeffF-coeffB) * ChiPrime * (FitC * FracMom2(FitV+1.0, FitS+2.0, moments) + 2.0 * moments[7] + 2.0 * FitC * FracMom2(FitV, FitS+2.0, moments) + moments[3] + FitC * FracMom2(FitV-1.0, FitS+2.0, moments)); //d(dc)=d(dp)
//		return (coeffF-coeffB) * ChiPrime * ((8.0/3.0) * moments[7] + (7.0/3.0) * moments[3] + (2.0/3.0) * FracMom2(-1.0, 2.0, moments)); //d(np)=0
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(8.0/3.0, moments) - FracMom2(2.0, 2.0, moments)); //Coalescence
	}
	else if (i==1 && j==2)
	{
//		return (coeffF-coeffB) * ChiPrime * (2.0 * FitC * FracMom2(FitV, FitS+3.0, moments) + FitC * FitC * FracMom2(2.0*FitV-1.0, 2.0*FitS+3.0, moments) + moments[6] + 2.0 * FitC * FracMom2(FitV-1.0, FitS+3.0, moments) + FitC * FitC * FracMom2(2.0*FitV-2.0, 2.0*FitS+3.0)); //d(dc)=d(dp)
//		return (coeffF-coeffB) * ChiPrime * ((7.0/3.0) * moments[6] + (16.0/9.0) * FracMom2(-1.0, 3.0, moments) + (4.0/9.0) * FracMom2(-2.0, 3.0, moments)); //d(np)=0
		return (coeffF-coeffB) * ChiPrime * (GetThetaSG(7.0/3.0, moments) - FracMom2(1.0, 3.0, moments)); //Coalescence
	}*/
	else
	{
		cerr << "#error: no need to compute source term for surface growth for M_" << i << "," << j << endl;
		exit( 2 );
	}

	return 0.0;
}

Double TSoot::SourceSurfGrowthNew( int i, Double *moments, Double *Y, Double density, Double *molarMass )
{
	Double	ChiPrime = 1.0e19 * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

	Double	fmom23 = FracMom2( 2.0 / 3.0, moments );

	Double	M0Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[0], 1.0e-60 ) );
	Double	M1Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[1], 1.0e-60 ) );

//	M0Fact = ( moments[0] < SMALLSOOT ) ? 0.0 : M0Fact;
//	M1Fact = ( moments[1] < SMALLSOOT ) ? 0.0 : M1Fact;

	Double	coeffF = GetSurfGrowthCoeffForw( Y, density, molarMass );
	Double	coeffB = M0Fact * M1Fact * GetSurfGrowthCoeffBackw( Y, density, molarMass ); //changed by GB
#ifdef OXMOM0NEW
	//Double	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 9.0 * SMALLSOOT ), 1.0 / 9.0 );changed by GB
//	Double	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 18.0 * SMALLSOOT ), 1.0 / 18.0 );
//	Double	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 18.0 * SMALLSOOT ), 1.0 );
#else
	Double	beta = fBeta;
#endif
//	Double	C_H = density * Y[f_H] / molarMass[f_H];


	switch( i ) {
		case 0: 
//			return 0.0;
		        return - coeffB * ChiPrime * FracMom2( -1.0 / 3.0, moments );
		case 1:
			return ( coeffF - coeffB ) * ChiPrime * fmom23;
		case 2:
			return ( coeffF - coeffB ) * ChiPrime * ( fmom23 + 2.0 * FracMom2( 5.0 / 3.0, moments ) );
		case 3:
			return ( coeffF - coeffB ) * ChiPrime * ( fmom23 + 3.0 * ( FracMom2( 5.0 / 3.0, moments ) 
						+ FracMom2( 8.0 / 3.0, moments ) ) );
		default:
			cerr << "#error: no need to compute source term for surface growth for M_" << i << NEWL;
			exit( 2 );

	}
	return 0.0;
}

Double TSoot::SourceSurfGrowthNew( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSurfGrowthCoeff( Y, density, molarMass );
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return coeff * fracMom[k02_3];
		case 2:
			return coeff * ( fracMom[k02_3] + 2.0 * fracMom[k05_3] );
		case 3:
			return coeff * ( fracMom[k02_3] + 3.0 * ( fracMom[k05_3] + fracMom[k08_3] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );

	}
	return 0.0;
}

Double TSoot::GetSurfGrowthCoeff( Double *Y, Double density, Double *molarMass )
{
	Double			C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
#ifdef NEWSURFGROWTH
	Double			C_H = density * Y[f_H] / molarMass[f_H];
	Double			*k = fSootRateCoeffs->vec;
	Double			locSG = k[ks10f] * fThirdBodyConc * C_C2H2 * fFk10 * fCSootStar 
							- k[ks13b] * C_H * ( 1.0 - fFk10 );

	return fAlpha * locSG;
#else
	return fAlpha * fCSootStar * fSootRateCoeffs->vec[ks10f] * C_C2H2;
#endif
}

Double TSoot::GetSurfGrowthCoeffForw( Double *Y, Double density, Double *molarMass )
{
	Double			C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
	Double			C_H = density * Y[f_H] / molarMass[f_H];
#ifdef NEWSURFGROWTH
	//Double			locSG = fSootRateCoeffs->vec[ks10f] * fThirdBodyConc * C_C2H2
	//				* fFk10 * fCSootStar;  changed by GB
	Double			locSG = (fSootRateCoeffs->vec[ks10f] * fThirdBodyConc * C_C2H2 
					* fCSootStar + fSootRateCoeffs->vec[ks13b] * C_H ) * fFk10;

//	return fAlpha * locSG;
	return locSG;
#else
	return fAlpha * fCSootStar * fSootRateCoeffs->vec[ks10f] * C_C2H2;
#endif
}

Double TSoot::GetSurfGrowthCoeffBackw( Double *Y, Double density, Double *molarMass )
{
#ifdef NEWSURFGROWTH
	Double			C_H = density * Y[f_H] / molarMass[f_H];
	//Double			locSG = fSootRateCoeffs->vec[ks13b] * C_H * ( 1.0 - fFk10 ); changed by GB
	Double			locSG = fSootRateCoeffs->vec[ks13b] * C_H;

//	return fAlpha * locSG;
	return locSG;
#else
	return 0.0;
#endif
}



//Mueller (8/20/07)
Double TSoot::GetThetaOx(Double frac, Double *moments)
{
	Double	theta0, theta1, theta2, theta3;

	theta0 = moments[1];
	theta1 = FracMom2(1.0, 1.0, moments) - moments[1];

	if (fNSootMoments > 3)
		theta2 = FracMom2(2.0, 1.0, moments) - 2.0 * moments[4] + moments[1];

	if (fNSootMoments > 6)
		theta3 = FracMom2(3.0, 1.0, moments) - 3.0 * moments[8] + 3.0 * moments[4] - moments[1];

	if (fNSootMoments==3)
		return pow(theta0, 1.0-frac) * pow(theta1*theta1, 0.5*frac);
	else if (fNSootMoments==6)
		return pow(theta0, (1.0-frac)*(1.0-0.5*frac)) * pow(theta1, frac*(2.0-frac)) * pow(theta2, 0.5*frac*(frac-1.0));
	else if (fNSootMoments==10)
		return pow(theta0, (1.0-frac)*(1.0-0.5*frac)*(1.0-(1.0/3.0)*frac)) * pow(theta1, frac*(2.0-frac)*(1.5-0.5*frac)) * pow(theta2, 0.5*frac*(frac-1.0)*(3.0-frac)) * pow(theta3, (1.0/3.0)*frac*(0.5*frac-0.5)*(frac-2.0));

	return 0.0;
}

//Mueller
Double TSoot::SourceSootOxidationNew(int i, int j, Double *moments, Double *Y, Double density, Double *molarMass)
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );

	Double	ChiPrime = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

	if (i==0 && j==0)
	{
//		return 0.0;
		return -coeff * ChiPrime * FracMom2(-1.0, 1.0, moments);
	}
	else if (i==1 && j==0)
	{
		return - coeff * ChiPrime * moments[1];
	}
	else if (i==0 && j==1)
	{
//		return -(2.0/3.0) * coeff * ChiPrime * FracMom2(-1.0, 2.0, moments); //d(np)=0
		return 0.0 * coeff * ChiPrime * (GetThetaOx(2.0/3.0, moments) - FracMom2(0.0, 2.0, moments)); //Coalescence
	}
	else if (i==2 && j==0)
	{
		return -coeff * ChiPrime * (2.0 * moments[4] - moments[1]);
	}
	else if (i==0 && j==2)
	{
//		return -coeff * ChiPrime * ((4.0/3.0) * FracMom2(-1.0, 3.0, moments) - (4.0/9.0) * FracMom2(-2.0, 3.0, moments)); //d(np)=0
		return coeff * ChiPrime * (GetThetaOx(4.0/3.0, moments) - FracMom2(0.0, 3.0, moments)); //Coalescence
	}
	else if (i==1 && j==1)
	{
//		return -coeff * ChiPrime * ((5.0/3.0) * moments[3] - (2.0/3.0) * FracMom2(-1.0, 2.0, moments)); //d(np)=0
		return coeff * ChiPrime * (GetThetaOx(5.0/3.0, moments) - FracMom2(1.0, 2.0, moments)); //Coalescence
	}
	else if (i==3 && j==0)
	{
		return -coeff * ChiPrime * (3.0 * moments[8] - 3.0 * moments[4] + moments[1]);
	}
	else if (i==0 && j==3)
	{
//		return -coeff * ChiPrime * (2.0 * FracMom2(-1.0, 4.0, moments) - (4.0/3.0) * FracMom2(-2.0, 4.0, moments) + (8.0/27.0) * FracMom2(-3.0, 4.0, moments));
		return coeff * ChiPrime * (moments[8] - 2.0 * moments[4] + moments[1] - FracMom2(0.0, 4.0, moments)); //Coalescence
	}
	else if (i==2 && j==1)
	{
//		return -coeff * ChiPrime * ((8.0/3.0) * moments[7] - (7.0/3.0) * moments[3] + (2.0/3.0) * FracMom2(-1.0, 2.0, moments)); //d(np)=0
		return coeff * ChiPrime * (GetThetaOx(8.0/3.0, moments) - FracMom2(2.0, 2.0, moments)); //Coalescence
	}
	else if (i==1 && j==2)
	{
//		return -coeff * ChiPrime * ((7.0/3.0) * moments[6] - (16.0/9.0) * FracMom2(-1.0, 3.0, moments) + (4.0/9.0) * FracMom2(-2.0, 3.0, moments));
		return coeff * ChiPrime * (GetThetaOx(7.0/3.0, moments) - FracMom2(1.0, 3.0, moments)); //Coalescence
	}
	else
	{
		cerr << "#error: no need to compute source term for oxidation for M_" << i << "," << j << endl;
		exit( 2 );
	}

	return 0.0;
}

Double TSoot::SourceSootOxidationNew( int i, Double *moments, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );	// F4
//	Double	beta;
	
/*	if ( moments[1] < SMALLSOOT || moments[0] < SMALLSOOT ) {*/
/*		coeff = 0.0;*/
/*	}*/
	
	Double	M0Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[0], 1.0e-60 ) );
	Double	M1Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[1], 1.0e-60 ) );

//	M0Fact = ( moments[0] < SMALLSOOT ) ? 0.0 : M0Fact;
//	M1Fact = ( moments[1] < SMALLSOOT ) ? 0.0 : M1Fact;

#ifdef OXMOM0NEW
	//beta = MIN( M0Fact * moments[0] / MAX( moments[1], 9.0 * SMALLSOOT ), 1.0 / 9.0 );changed by GB
//	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 18.0 * SMALLSOOT ), 1.0 / 18.0 );
#else
	beta = fBeta;
#endif


	switch( i ) {
		case 0: 
		  //return -M0Fact * M1Fact * beta * coeff * FracMom2( -1.0 / 3.0, moments ); changed by GB
			return -coeff * FracMom2( -1.0 / 3.0, moments );
		case 1:
		  //return -M0Fact * M1Fact * coeff * FracMom2( 2.0 / 3.0, moments ); changed by GB
			return - coeff * FracMom2( 2.0 / 3.0, moments );
		case 2:
			return -coeff * ( -FracMom2( 2.0 / 3.0, moments ) 
								+ 2.0 * FracMom2( 5.0 / 3.0, moments ) );
		case 3:
			return 0.0/*-coeff * ( FracMom2( 2.0 / 3.0, moments ) 
								+ 3.0 * ( -FracMom2( 5.0 / 3.0, moments ) 
								+ FracMom2( 8.0 / 3.0, moments ) ) )*/;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::SourceSootOxidationNew( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );	// F4
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
#ifdef OXMOM0NEW
			fprintf( stderr, "#warning: not yet implemented\n" );
			return -fracMom[km1_3] / fracMom[k02_3] * coeff * fracMom[k02_3];
#else
			return -fBeta * coeff * fracMom[km1_3];
#endif

		case 1:
			return -coeff * fracMom[k02_3];
		case 2:
			return -coeff * ( ( 1.0 - fBeta ) * fracMom[k02_3] + 2.0 * fracMom[k05_3] );
		case 3:
			return -coeff * ( ( 1.0 - fBeta ) * ( fracMom[k02_3] + 3.0 * fracMom[k05_3] ) 
								+ 3.0 * fracMom[k08_3] );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

/*Double TSoot::SourceSootOxidation( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );	// F4
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
			return -fBeta * coeff * fracMom[km1_24];
		case 1:
			return -coeff * fracMom[k23_24];
		case 2:
			return -coeff * ( ( 1.0 - fBeta ) * fracMom[k23_24] + 2.0 * fracMom[k47_24] );
		case 3:
			return -coeff * ( ( 1.0 - fBeta ) * ( fracMom[k23_24] + 3.0 * fracMom[k47_24] ) 
								+ 3.0 * fracMom[k71_24] );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );

	}
}
*/

Double TSoot::GetSootOxCoeff( Double *Y, Double density, Double *molarMass )
{
	Double			C_OH = density * Y[f_OH] / molarMass[f_OH];
	Double			C_O2 = density * Y[f_O2] / molarMass[f_O2];
	Double                  C_H = density * Y[f_H] / molarMass[f_H]; //added by GB
	Double                  C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2]; //added by GB
	Double                  fCSootC2H2;
	Double			*k = fSootRateCoeffs->vec;

#ifdef NEWSURFGROWTH

	fCSootC2H2 = ( k[ks10f] * C_C2H2 * fThirdBodyConc * fCSootStar + k[ks13b] * C_H ) * fFk10 / k[ks13f]; //added by GB
	//return fAlpha * ( k[ks111] * C_O2 * fCSootStar + k[ks12] * C_OH ); changed by GB
//	return fAlpha * ( ( k[ks111] * fCSootStar + k[ks112] * fCSootC2H2 ) * C_O2 + k[ks12] * C_OH );
	return ((k[ks111] * fCSootStar + k[ks112] * fCSootC2H2) * C_O2 + k[ks12] * C_OH);

#else
	return fAlpha * ( k[ks10b] * ( density * Y[f_H] / molarMass[f_H] ) 
				+ k[ks12] * C_OH + k[ks11] * C_O2 * fCSootStar );
#endif
}

/*
void T1DSoot::FillJacobi( T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate )
{
	int				i, ioff, l;
	int				tempOff = flame->GetOffsetTemperature();
	int				specOff = flame->GetOffsetFirstSpecies();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			*pahMoments = flameNode->pahMoments;
	Double			**a = nodeInfo->a;
	Double			hnenn = nodeInfo->hnenn;
	Double			temp = flameNode->temp[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;
	
	ComputeFractionalMoments( moments );
	ComputePhi( temp );

	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );

	if ( coordinate == kPhysical ) {
		int fVVelocity = flame->GetOffsetVVelocity();
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i+fOffsetMoments;

			//	convection
			FillJacNonlinearConvectUpwind( fVVelocity, ioff, nodeInfo, 1.0 );

			//	diffusion
			if ( fSizeDepDiff ) {
#ifndef NODIFF
				FillJacSootDiffusion( i, -1.0, kPhysical, flame, nodeInfo );
#endif
			}
			else {
				FillJacSootDiffusionNew( i, -1.0, kPhysical, flame, nodeInfo );
			}

			// thermophoresis
			if ( fThermoPhoresis ) {
				FillJacSootThermoPhoresis( i, -1.0, kPhysical, flame, nodeInfo );
			}

			// source term m_i/rho
			for ( l = 0; l < fNSootMoments; ++l ) {
				a[l+fOffsetMoments][ioff] -= flameNode->dMdx[l][i] * hnenn;
			}
			a[tempOff][ioff] -= flameNode->dMdx[fNSootMoments][i] * hnenn;
			int	nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
			for ( l = 0; l < nSpeciesIn; ++l ) {
				a[l+specOff][ioff] -= flameNode->dMdx[l+1+fNSootMoments][i] * hnenn;
			}

/*			// nucleation
			if ( fNucleation ) {
				a[tempOff][ioff] -= 0.5 / temp * NucleationNew( i, temp, pahMoments )
								 / density * hnenn;
			}

			// coagulation
			if ( fCoagulation ) {
				// 1.0 / density already in 'SootCoagFunc' (since y is used for moments)
				a[tempOff][ioff] -= dfdyUpwind( tempOff, ioff, SootCoagFunc, nodeInfo, flame )
											* hnenn;
//				a[ioff][ioff] -= FillJacSourceCoagulation( i, i, moments )
//											* hnenn;
//				a[ioff][ioff] -= dfdyUpwind( ioff, ioff, SootCoagFunc, nodeInfo, flame )
//											* hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] -= dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame )
												* hnenn;
//					a[l+fOffsetMoments][ioff] -= FillJacSourceCoagulation( i, l, moments )
//												* hnenn;
				}
			}

			//	condensation
			if ( fCondensation ) {
				a[tempOff][ioff] -= SourceCondensationNew( i, temp, pahMoments
									, moments, Y, density, molarMass ) * hnenn / ( 2.0 * temp * density );
				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] -= dfdyUpwind( l+fOffsetMoments, ioff, SootCondFunc, nodeInfo, flame )
//												* hnenn;
//					a[l+fOffsetMoments][ioff] -= FillJacCondensationNew( l, i, temp, pahMoments, moments )
//												/ density * hnenn;
				}
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
//				FillJacSurfGrowthNew( i, -1.0 / density, nodeInfo, flame );
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
//				FillJacSootOxidationNew( i, -1.0 / density, nodeInfo, flame );
//				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] -= dfdyUpwind( l+fOffsetMoments, ioff, SootOxidationFunc, nodeInfo, flame )
//												* hnenn;
//				}
			}*//*
		}
	}
	else if ( coordinate == kSimilarity ){
		int fVVelocity = flame->GetOffsetVVelocity();
		Double	oneOverRhoMuRef = 1.0 / ( flameNode->rhoInf 
						* flameNode->viscosityInf );
		Double	oneOverRhoA = 1.0 / ( flameNode->mixDensity[kCurr] 
						* flame->GetStrainRate() );
						
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i+fOffsetMoments;

			//	convection
			FillJacNonlinearConvectUpwind( fVVelocity, ioff, nodeInfo, 1.0, FALSE );

			//	diffusion
			if ( fSizeDepDiff ) {
				FillJacSootDiffusion( i, oneOverRhoMuRef, kSimilarity, flame, nodeInfo );
			}
			else {
				FillJacSootDiffusionNew( i, oneOverRhoMuRef, kSimilarity, flame, nodeInfo );
			}

			// thermophoresis
			if ( fThermoPhoresis ) {
				FillJacSootThermoPhoresis( i, oneOverRhoMuRef, kSimilarity, flame, nodeInfo );
			}

			// nucleation
			if ( fNucleation ) {
				a[tempOff][ioff] += 0.5 * oneOverRhoA * NucleationNew( i, temp, pahMoments ) 
												/ temp * hnenn;
			}

			// coagulation
			if ( fCoagulation ) {
				a[tempOff][ioff] += oneOverRhoA * dfdyUpwind( tempOff, ioff, SootCoagFunc, nodeInfo, flame )
										* hnenn;
//				a[ioff][ioff] += oneOverRhoA * FillJacSourceCoagulation( i, i, moments )
//											* hnenn;
//				a[ioff][ioff] += oneOverRhoA * dfdyUpwind( ioff, ioff, SootCoagFunc, nodeInfo, flame )
//											* hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += oneOverRhoA * dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame )
												* hnenn;
//					a[l+fOffsetMoments][ioff] += oneOverRhoA * FillJacSourceCoagulation( i, l, moments )
//												* hnenn;
#ifdef DEBUGCOAG
					Double	coagAnalyt = FillJacSourceCoagulation( i, l, moments );
					Double	coagNum = dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame );
					if ( fabs( ( coagAnalyt - coagNum ) / coagNum ) > 0.001 ) {
						fprintf( stderr, "gp = %d\tdCoag_%d/dM_%d Num = %g\tAnalyt= %g\n"
									, nodeInfo->gridPoint, i, l, coagNum, coagAnalyt );
					}
#endif
				}
			}

			//	condensation
			if ( fCondensation ) {
//				FillJacSourceCondensation( i, oneOverRhoA, nodeInfo, flame );
				a[tempOff][ioff] += oneOverRhoA * dfdyUpwind( tempOff, ioff, SootCondFunc, nodeInfo, flame )
											* hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] += oneOverRhoA * dfdyUpwind( l+fOffsetMoments, ioff, SootCondFunc, nodeInfo, flame )
//												* hnenn;
					a[l+fOffsetMoments][ioff] += density * oneOverRhoA * FillJacCondensationNew( l, i, temp, pahMoments, moments )
												* hnenn;
				}
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
				FillJacSurfGrowthNew( i, oneOverRhoA, nodeInfo, flame );
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
				FillJacSootOxidationNew( i, oneOverRhoA, nodeInfo, flame );
/*				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += oneOverRhoA * dfdyUpwind( l+fOffsetMoments, ioff, SootOxidationFunc, nodeInfo, flame )
												* hnenn;
				}*//*
			}
		}
	}
	else if ( coordinate == kMixtureFraction ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i+fOffsetMoments;

			//	diffusion
			if ( fSizeDepDiff ) {
//				FillJacSecondDerivCentral( ioff, ioff, 1.0, nodeInfo );
//				FillJacSootDiffusion( i, 1.0 / ( 2.0 * GetLewis1() ), kMixtureFraction, flame, nodeInfo );
			}
			else {
//				FillJacSootDiffusionNew( i, 1.0 / ( 2.0 * GetLewis1() ), kMixtureFraction, flame, nodeInfo );
			}

			// nucleation
			if ( fNucleation ) {
				a[tempOff][ioff] += 0.5 / density * NucleationNew( i, temp, pahMoments ) / temp * hnenn;
			}

			// coagulation
			if ( fCoagulation ) {
				a[tempOff][ioff] += dfdyUpwind( tempOff, ioff, SootCoagFunc, nodeInfo, flame )
											 / density * hnenn;
//				a[ioff][ioff] += FillJacSourceCoagulation( i, i, moments )
//											/ density * hnenn;
//				a[ioff][ioff] += dfdyUpwind( ioff, ioff, SootCoagFunc, nodeInfo, flame )
//											/ density * hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame )
												/ density * hnenn;
//					a[l+fOffsetMoments][ioff] += FillJacSourceCoagulation( i, l, moments )
//												/ density * hnenn;
				}
			}

			//	condensation
			if ( fCondensation ) {
				a[tempOff][ioff] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) 
						/ ( 2.0 * density * temp ) * hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += FillJacCondensationNew( l, i, temp, pahMoments, moments )
												 / density * hnenn;
				}
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
				FillJacSurfGrowthNew( i, 1.0 / density, nodeInfo, flame );
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
				FillJacSootOxidationNew( i, 1.0 / density, nodeInfo, flame );
//				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] += dfdyUpwind( l+fOffsetMoments, ioff, SootOxidationFunc, nodeInfo, flame )
//												/ density * hnenn;
//				}
			}
		}
	}
}*/

/*void TSoot::InitSootReactions( void )
{
	Double	*a = fASoot->vec;
	Double	*n = fNSoot->vec;
	Double	*eOverR = fEOverRSoot->vec;
	
	a[ks8f] = 7.900e10;
	n[ks8f] = 0.0;
	eOverR[ks8f] = 41.70e6 / RGAS;

	a[ks8b] = 5.683e9;
	n[ks8b] = 0.0;
	eOverR[ks8b] = 20.67e6 / RGAS;

	a[ks9] = 1.000e10;
	n[ks9] = 0.0;
	eOverR[ks9] = 0.0;

	a[ks10f] = 1.000e10;
	n[ks10f] = 0.0;
	eOverR[ks10f] = 0.0;

	a[ks10b] = 2.000e18;
	n[ks10b] = 0.0;
	eOverR[ks10b] = 290.9e6 / RGAS;

	if ( fSurfaceOxidation ) {
		a[ks11] = 1.0e10;
		n[ks11] = 0.0;
		eOverR[ks11] = 0.0;
	
		a[ks12] = 1.3e10;
		n[ks12] = 0.0;
		eOverR[ks12] = 46.0e6 / RGAS;
	}
	else {
		a[ks11] = 0.0;
		n[ks11] = 0.0;
		eOverR[ks11] = 0.0;
	
		a[ks12] = 0.0;
		n[ks12] = 0.0;
		eOverR[ks12] = 0.0;
	}
}
*/

//Additional coagulation source terms used for coalescence study.
/*
Double TSoot::SourceCoagulationNew(int i, int j, Double temp, Double *moments, Double nuc)
{
	Double coag = -0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments);

	if (i==0 && j==0)
	{
//		return -0.5 * GetPhi(0, 0, 0, 0, temp, 1.8, moments, moments); //Aggregation
		return -0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==1 && j==0)
	{
		return 0.0;
	}
	else if (i==0 && j==1)
	{
//		return nuc * ((1.0/3.0) * pow(moments[0], -2.0/3.0) * pow(moments[2], 2.0/3.0) + (2.0/3.0) * 18.0 * pow(moments[0], 1.0/3.0) * pow(moments[2], -1.0/3.0) - pow(18.0, 2.0/3.0)) + coag * (1.0/3.0) * pow(moments[0], -2.0/3.0) * pow(moments[2], 2.0/3.0); //Explicit source in volume
		return nuc * ((1.0/3.0) * pow(moments[0], -1.0/1.0) * pow(moments[1], 1.0/1.0) + (2.0/3.0) * 18.0 * pow(moments[0], 1.0/2.0) * pow(moments[1], -1.0/2.0) - pow(18.0, 2.0/3.0)) + coag * (1.0/3.0) * pow(moments[0], -1.0/1.0) * pow(moments[1], 1.0/1.0); //Explicit source in surface
	}
	else
	{
		cerr << "#error: no need to compute source term for coagulation for M_" << i << "," << j << endl;
		exit( 2 );
	}

	return 0.0;
}

//Mueller (7/20/07)
Double TSoot::SourceCoagulationNew(int i, int j, Double temp, Double *moments)
{
	if (i==0 && j==0)
	{
//		return -0.5 * GetPhi(0, 0, 0, 0, temp, 1.8, moments, moments); //Aggregation
		return -0.5 * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments); //Coalescence
	}
	else if (i==1 && j==0)
	{
		return 0.0;
	}
	else if (i==0 && j==1)
	{
//		return 0.0; //Aggregation
//		return (1.0/2.0) * pow(GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments), 1.0/3.0) * pow(2.0 * GetPhi(1, 0, 0, 0, temp, 3.0, moments, moments), 2.0/3.0) - GetPhi(2.0/3.0, 0, 0, 0, temp, 3.0, moments, moments); //Coalescence with only volume
//		return (1.0/2.0) * pow(GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments), 1.0/3.0) * pow(2.0 * GetPhi(1, 0, 0, 0, temp, 3.0, moments, moments), 2.0/3.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence with both surface and volume
//		return (1.0/2.0) * pow(GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments), 1.0/3.0) * pow(2.0 * GetPhi(0, 1.5, 0, 0, temp, 3.0, moments, moments), 2.0/3.0) - GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments); //Coalescence with only surface
//		return (1.0/2.0) * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments, (7.0/6.0)) - GetPhi(2.0/3.0, 0, 0, 0,  temp, 3.0, moments, moments); //Better Interpolation with volume
		return (1.0/2.0) * GetPhi(0, 0, 0, 0, temp, 3.0, moments, moments, (7.0/6.0)) - GetPhi(0, 1, 0, 0,  temp, 3.0, moments, moments); //Better Interpolation with volume and surface
//		return ((1.0/2.0) * pow(2.0, 2.0/3.0) - 1.0) * GetPhi(2.0/3.0, 0, 0, 0, temp, 3.0, moments, moments); //Fewer terms with volume
//		return ((1.0/2.0) * pow(2.0, 2.0/3.0) - 1.0) * GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments); //Fewer terms with surface
//		return GetC(temp, 3.0) * ((1.0/2.0) * pow(2.0, 2.0/3.0) - 1.0) * pow(2.0, 5.0/2.0) * pow(moments[2]/moments[0], 5.0/6.0) * moments[0] * moments[0]; //Assumption of delta function with volume
//		return GetC(temp, 3.0) * ((1.0/2.0) * pow(2.0, 2.0/3.0) - 1.0) * pow(2.0, 5.0/2.0) * pow(moments[1]/moments[0], 5.0/4.0) * moments[0] * moments[0]; //Assumption of delta function with surface
//		return (1.0/3.0) * GetPhi(-1, 1, 1, 0, temp, 3.0, moments, moments) - (1.0/2.0) * GetPhi(0, 1, 0, 0, temp, 3.0, moments, moments); //dS/S=(2/3)*dV/V
	}
	else
	{
		cerr << "#error: no need to compute source term for coagulation for M_" << i << "," << j << endl;
		exit( 2 );
	}

	return 0.0;
}
*/


/*Double SootCoagFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag )
{
	T1DFlamePtr 	flame = ( T1DFlamePtr ) object;
//	TFlameNodePtr	flameNode = flame->GetFlameNode();
	T1DSootPtr		soot = flame->GetSoot();
	int				momOff = soot->GetOffsetSootMoments();
	Double			*phi = soot->fPhi->vec;
	Double			temp = nodeInfo->y[flame->GetOffsetTemperature()];
//	Double			temp = flameNode->temp[kCurr];
	Double			*moments = &nodeInfo->y[momOff];
	
	soot->ComputeFractionalMoments( moments );
	soot->ComputePhi( temp );
	
	switch ( j - momOff ) {
		case 0:
			return -0.5 * phi[kh00];
		case 1:
			return 0.0;
		case 2:
			return phi[kh11];
		case 3:
			return 3.0 * phi[kh12];
		default:
			cerr << "#error: no need to compute coagulation source_" << j - momOff << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double SootCondFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag )
{
	T1DFlamePtr 	flame = ( T1DFlamePtr ) object;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	T1DSootPtr		soot = flame->GetSoot();
	int				momOff = soot->GetOffsetSootMoments();
	Double			density = flameNode->mixDensity[kCurr];
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			temp = nodeInfo->y[flame->GetOffsetTemperature()];
//	Double			temp = flameNode->temp[kCurr];
	Double			*mom = &nodeInfo->y[momOff];
	Double			*pahMom = flame->GetFlameNode()->pahMoments;
		
	return soot->SourceCondensationNew( j - momOff, temp, pahMom, mom, Y, density, molarMass );
}

Double SootOxidationFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag )
{
	T1DFlamePtr 	flame = ( T1DFlamePtr ) object;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			source;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	T1DSootPtr		soot = flame->GetSoot();
	Double			temp = flameNode->temp[kCurr];
	int				momOff = soot->GetOffsetSootMoments();
	Double			*kSoot = soot->fSootRateCoeffs->vec;
	int				nSootMoments = soot->GetNSootMoments();
	Double			*moments = New1DArray( nSootMoments );

	for ( int i = 0; i < nSootMoments; ++i ) {
		moments[i] *= density;
	}

	soot->ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	source = soot->SourceSootOxidationNew( j-momOff, moments, Y, density, molarMass );
	Free1DArray( moments );
	return source;
}*/
/*
void T1DSoot::PrintSurfDepCoag( TNewtonPtr bt, T1DFlamePtr flame )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	FILE				*fp = NULL;
	Double		dummy;
	int			counter;
	
	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sSurfDepCoag%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) { 
		cerr << "#warning: unable to open file " << flame->GetOutFileBuff() << NEWL;
		exit(2);
	}
	
	fprintf( fp, "*\n%-12s", "y" );
	fprintf( fp, "\tratio\tCoagNew/Coag\tksg\tc*phi\tc\tS\ttemp\tcoag\tCst\n" );
		
	for ( k = 0; k < N; ++k ){
		bt->SetNodeInfo( flame, k );
		PrintSurfDepCoag( flame, nodeInfo, fp );
	}
    fclose( fp );
}

void T1DSoot::PrintSurfDepCoag( T1DFlamePtr flame, NodeInfoPtr nodeInfo, FILE *fp )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			temp = flameNode->temp[kCurr];
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
	Double	ksg = GetSurfGrowthCoeff( Y, density, molarMass );
	Double	c = GetC( temp );
	Double	coag = -0.5 * fPhi->vec[kh00];
	Double	S = GetSCoag( moments, temp );
	Double	Cst = 1.0 / ( 6.0 * fCoagFact );
	Double	newCoag = coag * ( 1.0 - 1.0 / ( 1.0 + Cst * fPhi->vec[kh00] * ksg / ( c * c * S ) ) );
	fprintf( fp, "%-.6e", *nodeInfo->x );


	fprintf( fp, "\t%-.6e", fPhi->vec[kh00]*Cst*ksg/(c*c*S) );
	fprintf( fp, "\t%-.6e", newCoag/coag );
	fprintf( fp, "\t%-.6e", ksg );
	fprintf( fp, "\t%-.6e", fPhi->vec[kh00] );
	fprintf( fp, "\t%-.6e", c );
	fprintf( fp, "\t%-.6e", S );
	fprintf( fp, "\t%-.6e", temp);
	fprintf( fp, "\t%-.6e", coag);
	fprintf( fp, "\t%-.6e", Cst );
	fprintf( fp, "\n" );
}*/


/*
Double T1DSoot::FillJacCondensationNew( int l, int i, Double temp, Double *pahMom, Double *mom )
{
	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return GetDerivPhiCond( l, 0, 1, temp, mom, pahMom );
		case 2:
			return GetDerivPhiCond( l, 0, 2, temp, mom, pahMom ) 
				+ 2.0 * GetDerivPhiCond( l, 1, 1, temp, mom, pahMom );
		case 3:
			return GetDerivPhiCond( l, 0, 3, temp, mom, pahMom ) 
				+ 3.0 * ( GetDerivPhiCond( l, 1, 2, temp, mom, pahMom ) 
						+ GetDerivPhiCond( l, 2, 1, temp, mom, pahMom ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

void T1DSoot::FillJacSurfGrowth( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int		iOff = i + fOffsetMoments;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	density = flameNode->mixDensity[kCurr];
	Double	**a = nodeInfo->a;
	Double	*moments = flameNode->moments;
	Double	*Y = flameNode->Y[kCurr];

	constCoeff *= GetSurfGrowthCoeff( Y, density, molarMass ) * nodeInfo->hnenn;

	switch( i ) {
		case 0: 
			break;
		case 1:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k23_24, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k23_24, moments );
			break;
		case 2:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k23_24, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k23_24, moments );

			a[fOffsetMoments+1][iOff] += 2.0 * constCoeff * GetDerivFracMom( 1, k47_24, moments );
			a[fOffsetMoments+2][iOff] += 2.0 * constCoeff * GetDerivFracMom( 2, k47_24, moments );
			break;
		case 3:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k23_24, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k23_24, moments );

			a[fOffsetMoments+1][iOff] += 3.0 * constCoeff * GetDerivFracMom( 1, k47_24, moments );
			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k47_24, moments );

			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k71_24, moments );
			a[fOffsetMoments+3][iOff] += 3.0 * constCoeff * GetDerivFracMom( 3, k71_24, moments );
			break;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
}

void T1DSoot::FillJacSurfGrowthNew( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int		iOff = i + fOffsetMoments;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	density = flameNode->mixDensity[kCurr];
	Double	**a = nodeInfo->a;
	Double	*moments = flameNode->moments;
	Double	*Y = flameNode->Y[kCurr];

	constCoeff *= GetSurfGrowthCoeff( Y, density, molarMass ) * nodeInfo->hnenn;

	switch( i ) {
		case 0: 
			break;
		case 1:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, 2.0/3.0, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, 2.0/3.0, moments );
//			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k02_3, moments );
//			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k02_3, moments );
			break;
		case 2:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] += 2.0 * constCoeff * GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] += 2.0 * constCoeff * GetDerivFracMom( 2, k05_3, moments );
			break;
		case 3:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] += 3.0 * constCoeff * GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k05_3, moments );

			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k08_3, moments );
			a[fOffsetMoments+3][iOff] += 3.0 * constCoeff * GetDerivFracMom( 3, k08_3, moments );
			break;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
}

void T1DSoot::FillJacSootRadiation( Double fact, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	// alphas taken from:
	//		Hubbard, G. L. , Tien C. L: 
	//		Infrared Mean Absorption Coefficient
	//		of Luminous Flames and Smoke
	//		Journal of Heat Transfer, vol 100, p. 235ff, 1978

	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int				tempOff = flame->GetOffsetTemperature();
	Double			temp = flameNode->temp[kCurr];
	const Double	alpha = -3.75e5, beta = 1735.0;
	Double			*moments = flameNode->moments;
	Double			**a = nodeInfo->a;
	Double			alphas = alpha + beta * temp;			// [m^-1] 
	Double			rad = GetSootRadiation( temp, moments );
	
	fact *= rad * nodeInfo->hnenn;

	a[fOffsetMoments+1][tempOff] += fact * flameNode->mixDensity[kCurr] / moments[1];

	a[tempOff][tempOff] += fact * ( 4.0 / temp 
			+ beta / alphas );
}

void T1DSoot::FillJacSootOxidationNew( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int		iOff = i + fOffsetMoments;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	density = flameNode->mixDensity[kCurr];
	Double	**a = nodeInfo->a;
	Double	*moments = flameNode->moments;
	Double	*Y = flameNode->Y[kCurr];
	Double	beta;

	constCoeff *= GetSootOxCoeff( Y, density, molarMass ) * nodeInfo->hnenn;

	switch( i ) {
		case 0: 
#ifdef OXMOM0NEW
			beta = MIN( ( ( moments[0] > SMALLSOOT ) ? moments[0] : 0.0 ) 
						/ MAX( moments[1], SMALLSOOT ), 1.0 / 9.0 );
#else
			beta = fBeta;
#endif
			a[fOffsetMoments][iOff] -= constCoeff * beta * GetDerivFracMom( 0, km1_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * beta * GetDerivFracMom( 1, km1_3, moments );
			break;
		case 1:
			a[fOffsetMoments][iOff] -= constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * GetDerivFracMom( 1, k02_3, moments );
			break;
		case 2:
			a[fOffsetMoments][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] -= 2.0 * constCoeff * GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] -= 2.0 * constCoeff * GetDerivFracMom( 2, k05_3, moments );
			break;
		case 3:
			a[fOffsetMoments][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] -= 3.0 * constCoeff * ( 1.0 - fBeta )
											* GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] -= 3.0 * constCoeff * ( 1.0 - fBeta )
											* GetDerivFracMom( 2, k05_3, moments );

			a[fOffsetMoments+2][iOff] -= 3.0 * constCoeff * GetDerivFracMom( 2, k08_3, moments );
			a[fOffsetMoments+3][iOff] -= 3.0 * constCoeff * GetDerivFracMom( 3, k08_3, moments );
			break;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
}*/

/*Double T1DSoot::FillJacSourceCoagulation( int i, int l, Double *moments )
{
//	returns dS_i/dM_l
		
	switch ( i ) {
		case 0:
			return -0.5 * GetDerivPhi( l, kh00, moments );
		case 1:
			return 0.0;
		case 2:
			return GetDerivPhi( l, kh11, moments );
		case 3:
			return 3.0 * GetDerivPhi( l, kh12, moments );
		default:
			cerr << "#error: no need to compute derivative of coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}*/

/*
void T1DSoot::FillJacSootThermoPhoresis( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
// fills the jacobian with     constCoeff * d/dy ( 0.55 * nu / T * M_r * dT/dy )

	const Double	A = 0.9;	// accommodation coefficient following
								// R. J. Santoro: The Transport and Growth
								// of Soot Particles ...
								// Comb. Sci. Tech., 1987, Vol. 53 p. 89
	static Double	pi = 4.0 * atan( 1.0 );
	const Double	fact = 3.0 / ( 4.0 * ( 1.0 + pi * A / 8.0 ) );
	constCoeff *= fact;
	int		rOff = r + fOffsetMoments;
	int		tempOff = flame->GetOffsetTemperature();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double	mu = flameNode->mixViscosity[kCurr];
	Double	muNext = flameNode->mixViscosity[kNext];
	Double	muPrev = flameNode->mixViscosity[kPrev];
	Double	M = nodeInfo->y[rOff];
	Double	MNext = nodeInfo->yNext[rOff];
	Double	MPrev = nodeInfo->yPrev[rOff];
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*temp = flameNode->temp;
	Double	*rho = flameNode->mixDensity;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	diffPlusHm, diffMinusH;
	Double	coeffPlus, coeffCurr, coeffMinus;

	if ( coordinate == kPhysical ) {
		diffPlusHm = hm * constCoeff 
				* ( mu * M / temp[kCurr]
					+ muNext * MNext / temp[kNext] );
		diffMinusH = h * constCoeff 
				* ( mu * M / temp[kCurr]
					+ muPrev * MPrev / temp[kPrev] );

		coeffPlus = constCoeff * muNext / temp[kNext]
					* hm * ( temp[kNext] - temp[kCurr] );
		coeffMinus = constCoeff * muPrev / temp[kPrev]
					* h * ( temp[kPrev] - temp[kCurr] );
		coeffCurr = constCoeff * mu / temp[kCurr]
					* ( hm * ( temp[kNext] - temp[kCurr] ) 
						+ h * ( temp[kPrev] - temp[kCurr] ) );
	}
	else {
		diffPlusHm = hm * constCoeff 
				* ( rho[kCurr] * mu * M / temp[kCurr]
					+ rho[kNext] * muNext * MNext / temp[kNext] );
		diffMinusH = h * constCoeff 
				* ( rho[kCurr] * mu * M / temp[kCurr]
					+ rho[kPrev] * muPrev * MPrev /temp[kPrev] );

		coeffPlus = constCoeff * rho[kNext] * muNext / temp[kNext]
					* hm * ( temp[kNext] - temp[kCurr] );
		coeffMinus = constCoeff * rho[kPrev] * muPrev / temp[kPrev]
					* h * ( temp[kPrev] - temp[kCurr] );
		coeffCurr = constCoeff * rho[kCurr] * mu / temp[kCurr]
					* ( hm * ( temp[kNext] - temp[kCurr] ) 
						+ h * ( temp[kPrev] - temp[kCurr] ) );
	}

	a[tempOff][rOff] -= diffPlusHm + diffMinusH;
	a[rOff][rOff] += coeffCurr;
	if ( !nodeInfo->lastPoint ) {
		b[tempOff][rOff] += diffPlusHm;
		b[rOff][rOff] += coeffPlus;
	}
	if ( !nodeInfo->firstPoint ) {
		c[tempOff][rOff] += diffMinusH;
		c[rOff][rOff] += coeffMinus;
	}
}*/

/*
void T1DSoot::PrintFracMom( T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*fm = fFracMom->vec;
	Double		**mo = fMoments->mat;
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sFracMoments_%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	fprintf( fp, "\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s"
				, "km1_3", "0", "k02_3", "1", "k05_3", "2", "k08_3", "3" );
	
	for ( k = 0; k < nGridPoints; ++k ) {
		ComputeFractionalMoments( mo[k] );
		fprintf( fp, "\n%-9E", x[k] );
		fprintf( fp, "\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E"
		, fm[km1_3], mo[k][0], fm[k02_3], mo[k][1], fm[k05_3], mo[k][2], fm[k08_3], mo[k][3] );
	}

	fclose( fp );
}

void T1DSoot::PrintPhi( T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*phi = fPhi->vec;
	Double		**mo = fMoments->mat;
	Double		*temp = flame->GetTemperature()->vec;
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sPhi_%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	fprintf( fp, "\tM_%-6s\tPhi_%-6s\tM_%-6s\tPhi_%-6s\tM_%-6s\tPhi_%-6s\tM_%-6s"
				, "0", "kh00", "1", "kh11", "2", "kh12", "3" );
	
	for ( k = 0; k < nGridPoints; ++k ) {
		ComputeFractionalMoments( mo[k] );
		ComputePhi( temp[k] );
		fprintf( fp, "\n%-9E", x[k] );
		fprintf( fp, "\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E"
		, mo[k][0], phi[kh00], mo[k][1], phi[kh11], mo[k][2], phi[kh12], mo[k][3] );
	}

	fclose( fp );
}

void T1DSoot::PrintDiffusivity( T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*diff = fSootDiff->vec;
	Double		**diffA1 = flame->GetSpecies()->GetDiffusivity()->mat;
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sdiff_%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	fprintf( fp, "\t%-6s\t%-6s"
				, "D_Soot", "D_A1" );
	
	fprintf( fp, "\n%-9E\t%-9E\t%-9E"
	, left, diff[kPrev], diffA1[kPrev][f_A1] );

	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		fprintf( fp, "\t%-9E\t%-9E"
		, diff[k], diffA1[k][f_A1] );
	}

	fprintf( fp, "\n%-9E\t%-9E\t%-9E"
	, right, diff[nGridPoints], diffA1[nGridPoints][f_A1] );

	fclose( fp );
}*/

/*void T1DSoot::FillJacSootConvection( int i, T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag velocityPositive )
{
//	velocityPositive should have the value FALSE, if 'V' has the negative 
//	direction of the physical velocity

// fills the jacobian with     V * d(M_i/rho)/dy

	
	int		iOff = i + fOffsetMoments;
	int		vOff = flame->GetOffsetVVelocity();
	int		tempOff = flame->GetOffsetTemperature();
	Double	*rho = flame->GetFlameNode()->mixDensity;
	Double	coeff;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;
	Double	**a = nodeInfo->a;
	Double	V = y[vOff];

	if ( ( V > 0.0 && velocityPositive ) || ( V < 0.0 && !velocityPositive ) ) {
		coeff = h * ( h + hm );
		a[vOff][iOff] += coeff * ( y[iOff] / rho[kCurr] - yPrev[iOff] / rho[kPrev] );
		a[iOff][iOff] += coeff * V / rho[kCurr];
//		a[tempOff][iOff] += coeff * V * y[iOff] / ( rho[kCurr] * y[tempOff] );
		if ( !nodeInfo->firstPoint ) {
			nodeInfo->c[iOff][iOff] -= coeff * V / rho[kPrev];
//			nodeInfo->c[tempOff][iOff] -= coeff * V * yPrev[iOff] / ( rho[kPrev] * yPrev[tempOff] );
		}
#ifdef FLUXBC
		else if ( !fSizeDepDiff ) {
			Double	c = rho[kPrev] * flame->GetFlameNode()->diffSoot[kPrev] / ( yPrev[vOff] * hm );
			a[iOff][iOff] -= V * coeff / rho[kCurr] * c / ( 1.0 + c ); 
		}
#endif
	}
	else {
		coeff = hm * ( h + hm );
		a[vOff][iOff] += coeff * ( yNext[iOff] / rho[kNext] - y[iOff] / rho[kCurr] );
		a[iOff][iOff] -= coeff * V / rho[kCurr];
//		a[tempOff][iOff] -= coeff * V * y[iOff] / ( rho[kCurr] * y[tempOff] );
		if ( !nodeInfo->lastPoint ) {
			nodeInfo->b[iOff][iOff] += coeff * V / rho[kNext];
//			nodeInfo->b[tempOff][iOff] += coeff * V * yNext[iOff] / ( rho[kNext] * yNext[tempOff] );
		}
#ifdef FLUXBC
		else if ( !fSizeDepDiff ) {
			Double	c = rho[kNext] * flame->GetFlameNode()->diffSoot[kNext] / ( yNext[vOff] * h );
			a[iOff][iOff] -= V * coeff / rho[kCurr] * c / ( 1.0 - c ); 
		}
#endif
	}
}
*/

/*
void T1DSoot::FillJacSootDiffusion( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
// fills the jacobian with     constCoeff * d/dy ( rho * diffusivity * d(M_(r-2/3)/rho)/dy )

	int		i;
	int		rOff = r + fOffsetMoments;
	int		tempOff = flame->GetOffsetTemperature();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
//	Double	diff = flameNode->diffusivity[f_A1];
//	Double	diffNext = flameNode->diffusivityNext[f_A1];
//	Double	diffPrev = flameNode->diffusivityPrev[f_A1];
	Double	diff = flameNode->diffSoot[kCurr];
	Double	diffNext = flameNode->diffSoot[kNext];
	Double	diffPrev = flameNode->diffSoot[kPrev];
	Double	*rho = flameNode->mixDensity;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	*moments = &y[fOffsetMoments];
	Double	*momentsNext = &yNext[fOffsetMoments];
	Double	*momentsPrev = &yPrev[fOffsetMoments];
	Double	fracIndex = r - 2.0 / 3.0;
	Double	fracMom = FracMom2( fracIndex, moments );
	Double	fracMomPrev = FracMom2( fracIndex, momentsPrev );
	Double	fracMomNext = FracMom2( fracIndex, momentsNext );
	Double	diffPlusHm, diffMinusH, fact, factNext, factPrev;

	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ){
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * constCoeff * nodeInfo->hm;
		diffMinusH = 2.0 * constCoeff * nodeInfo->h;
	}

	fact = ( diffPlusHm + diffMinusH ) * fracMom;
	factNext = diffPlusHm * fracMomNext;
	factPrev = diffMinusH * fracMomPrev;
	
#define FULLSOOTDIFFJAC

#ifdef FULLSOOTDIFFJAC
	for ( i = 0; i < fNSootMoments; ++i ) {
		a[i+fOffsetMoments][rOff] -= fact * GetAlphaI2( i, fracIndex ) 
											* rho[kCurr] / moments[i];
		if ( !nodeInfo->lastPoint ) {
			b[i+fOffsetMoments][rOff] += factNext * GetAlphaI2( i, fracIndex ) 
											* rho[kNext] / momentsNext[i];
		}
		if ( !nodeInfo->firstPoint ) {
			c[i+fOffsetMoments][rOff] += factPrev * GetAlphaI2( i, fracIndex ) 
											* rho[kPrev] / momentsPrev[i];
		}
	}
#else
	a[rOff][rOff] -= fact * GetAlphaI2( r, fracIndex ) * rho[kCurr] / moments[r];
	if ( !nodeInfo->lastPoint ) {
		b[rOff][rOff] += factNext * GetAlphaI2( r, fracIndex ) * rho[kNext] / momentsNext[r];
	}
	if ( !nodeInfo->firstPoint ) {
		c[rOff][rOff] += factPrev * GetAlphaI2( r, fracIndex ) * rho[kPrev] / momentsPrev[r];
	}
#endif
}

void T1DSoot::FillJacSootDiffusionNew( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
// fills the jacobian with     constCoeff * d/dy ( rho * diffusivity * d(M_r/rho)/dy )

	int		rOff = r + fOffsetMoments;
	int		tempOff = flame->GetOffsetTemperature();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double	diff = flameNode->diffSoot[kCurr];
	Double	diffNext = flameNode->diffSoot[kNext];
	Double	diffPrev = flameNode->diffSoot[kPrev];
	Double	*rho = flameNode->mixDensity;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	diffPlusHm, diffMinusH, fact, factNext, factPrev;

	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ){
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * constCoeff * nodeInfo->hm;
		diffMinusH = 2.0 * constCoeff * nodeInfo->h;
	}

	fact = ( diffPlusHm + diffMinusH );
	factNext = diffPlusHm;
	factPrev = diffMinusH;
	
	a[rOff][rOff] -= fact;
	if ( !nodeInfo->lastPoint ) {
		b[rOff][rOff] += factNext;
	}
#ifdef FLUXBC
	else if ( !fSizeDepDiff ) {
		Double	c = rho[kNext] * diffNext / ( yNext[flame->GetOffsetVVelocity()] * nodeInfo->h );
		a[rOff][rOff] -= factNext * rho[kPrev]/rho[kCurr] * c / ( 1.0 - c ); 
	}
#endif
	if ( !nodeInfo->firstPoint ) {
		c[rOff][rOff] += factPrev;
	}
#ifdef FLUXBC
	else if ( !fSizeDepDiff ) {
		Double	c = rho[kPrev] * diffPrev / ( yPrev[flame->GetOffsetVVelocity()] * nodeInfo->hm );
		a[rOff][rOff] += factPrev * rho[kPrev]/rho[kCurr] * c / ( 1.0 + c ); 
	}
#endif
}*/

/*void T1DSoot::FillRHS( T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate )
{

        int			i, ioff;
	TFlameNodePtr		flameNode = flame->GetFlameNode();
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			temp = flameNode->temp[kCurr];
	Double			density = flameNode->mixDensity[kCurr];
	Double			*rhs = nodeInfo->rhs;
	Double			*pahMoments = flameNode->pahMoments;
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;
	Double			*theSource = New1DArray( fNSootMoments );

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	
	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
	
	// convective term 
	if ( coordinate == kPhysical ) {
		int fVVelocity = flame->GetOffsetVVelocity();
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;

			//	convection
			rhs[ioff] += NonlinearConvectUpwind( nodeInfo->y[fVVelocity]
					, nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
					, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h );

			//	diffusion
			if ( fSizeDepDiff ) {
#ifndef NODIFF
				rhs[ioff] -= SootDiffusion( i, kPhysical, flame, nodeInfo );
#endif
			}
			else {
				rhs[ioff] -= SootDiffusionNew( i, kPhysical, flame, nodeInfo );
			}

//			Double diff = SootDiffusion( i, kPhysical, flame, nodeInfo );
//			Double tp = SootThermoPhoresis( i, kPhysical, flame, nodeInfo );
//			if ( fabs( diff ) * 10.0 < fabs( tp ) ) {
//				cerr << "T = " << temp << TAB
//					<< "Diff = " << diff
//					<< TAB << "ThermoPh = " << tp << NEWL;
//			}
//			if ( fabs( diff ) > fabs( tp ) * 10.0 ) {
//				cerr << TAB << "T = " << temp << TAB
//					<< "Diff = " << diff
//					<< TAB << "ThermoPh = " << tp << NEWL;
//			}

			// thermophoresis
			if ( fThermoPhoresis ) {
				rhs[ioff] -= SootThermoPhoresis( i, kPhysical, flame, nodeInfo );
			}
			rhs[ioff] -= flameNode->sootSource[i];

			// nucleation
//			if ( fNucleation ) {
//				theSource[i] += NucleationNew( i, temp, pahMoments ) / density;
//				rhs[ioff] -= NucleationNew( i, temp, pahMoments ) / density;
//			rhs[ioff] -= 0.5 * Y[f_A3R5AC] / molarMass[f_A3R5AC];
//			}

			//	coagulation
//			if ( fCoagulation ) {
//				theSource[i] += SourceCoagulation( i ) / density;
//				rhs[ioff] -= SourceCoagulation( i ) / density;
//			}

//			//	condensation
//			if ( fCondensation ) {
//				theSource[i] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;
//				rhs[ioff] -= SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;
//			}

			//	surface growth
//			if ( fSurfaceGrowth ) {
//				theSource[i] += SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;
//				rhs[ioff] -= SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;
//			}

			//	surface goxidation
//			if ( fSurfaceOxidation ) {
//				theSource[i] += SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;
//				rhs[ioff] -= SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;
//			}
//			if ( theSource[i] > 1.0 && fabs( theSource[i] - flameNode->sootSource[i] ) / MAX( fabs( theSource[i] ), 1.0e-30 ) > 1.0e-8 ) {
//				fprintf( stderr, "gp = %d\ti = %d\ttheSource = %g\tsootSource = %g\n"
//					, nodeInfo->gridPoint, i, theSource[i], flameNode->sootSource[i] );
//			}
		}
	}
	else if ( coordinate == kSimilarity ) {
		int fVVelocity = flame->GetOffsetVVelocity();
		Double	oneOverRhoMuRef = 1.0 / ( flameNode->rhoInf 
								* flameNode->viscosityInf );
		Double	oneOverRhoA = 1.0 / ( flameNode->mixDensity[kCurr] 
								* flame->GetStrainRate() );
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;
			//	convection
			rhs[ioff] += NonlinearConvectUpwind( nodeInfo->y[fVVelocity]
					, nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
					, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h, FALSE );

			//	diffusion
			if ( fSizeDepDiff ) {
				rhs[ioff] += oneOverRhoMuRef * SootDiffusion( i, kSimilarity, flame, nodeInfo );
			}
			else {
				rhs[ioff] += oneOverRhoMuRef * SootDiffusionNew( i, kSimilarity, flame, nodeInfo );
			}
			// thermophoresis
			if ( fThermoPhoresis ) {
				rhs[ioff] += oneOverRhoMuRef * SootThermoPhoresis( i, kSimilarity, flame, nodeInfo );
			}

			// nucleation
			if ( fNucleation ) {
				rhs[ioff] += oneOverRhoA * NucleationNew( i, temp, pahMoments );
			}
			//	coagulation
			if ( fCoagulation ) {
				fprintf( stderr, "###Attention: old Coagulation used here\n" );
				rhs[ioff] += oneOverRhoA * SourceCoagulation( i );
			}
			//	condensation
			if ( fCondensation ) {
				rhs[ioff] += oneOverRhoA * SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass );
			}
			//	surface growth
			if ( fSurfaceGrowth ) {
				rhs[ioff] += oneOverRhoA * SourceSurfGrowthNew( i, Y, density, molarMass );
			}
			//	surface oxidation
			if ( fSurfaceOxidation ) {
				rhs[ioff] += oneOverRhoA * SourceSootOxidationNew( i, moments, Y, density, molarMass );
			}
		}
	}
	else if ( coordinate == kMixtureFraction ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;

			//	diffusion
			if ( fSizeDepDiff ) {
//				rhs[ioff] += SecondDeriv( nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
//								, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h );
//				rhs[ioff] += SootDiffusion( i, kMixtureFraction, flame, nodeInfo )
//							/ ( 2.0 * GetLewis1() );
			}
			else {
//				rhs[ioff] += SootDiffusionNew( i, kMixtureFraction, flame, nodeInfo )
//						/ ( 2.0 * GetLewis1() );
			}

			// nucleation
			if ( fNucleation ) {
//				fprintf( stderr, "Nuc = %g\n", NucleationNew( i, temp, pahMoments ) );
				rhs[ioff] += NucleationNew( i, temp, pahMoments ) / density;
			}

			//	coagulation
			if ( fCoagulation ) {
//				rhs[ioff] += SourceCoagulation( i ) / density;
				rhs[ioff] += SourceCoagulationNew( i, temp, moments ) / density;
			}

			//	condensation
			if ( fCondensation ) {
				rhs[ioff] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
//				rhs[ioff] += SourceSurfGrowthNew( i, Y, density, molarMass ) / density;
				rhs[ioff] += SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;
			}

			//	surface goxidation
			if ( fSurfaceOxidation ) {
				rhs[ioff] += SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;
			}
		}
	}
	else {
		fprintf( stderr, "#error: invalid coordinate type %d\n", coordinate );
		exit( 2 );
	}
	Free1DArray( theSource );
}*/


// double TSoot::CoagQuad(double x, double y, double * moments, double temp)
// {
//   double C = GetC(temp);

//   //Compute elements of covariance matrix

//   double ux = log(moments[1]*moments[1]/pow(moments[0]*moments[0]*moments[0]*moments[3],1.0/2.0));
//   double uy = log(moments[2]*moments[2]/pow(moments[0]*moments[0]*moments[0]*moments[5],1.0/2.0));
//   double sx = sqrt(MAX(log(moments[3]*moments[0]/(moments[1]*moments[1])),0.0));
//   double sy = sqrt(MAX(log(moments[5]*moments[0]/(moments[2]*moments[2])),0.0));
//   double rho = log(moments[4]*moments[0]/(moments[1]*moments[2])) / (sx*sy);

//   //Check for delta function distribution
//   if (sx == 0 || sy == 0)
//     return C * pow(2.0,2.5) * pow(moments[1]/moments[0],1.0/6.0) * ((1.0/2.0) * pow(2.0*moments[1]/moments[0],x) * pow(2.0*moments[2]/moments[0],y) - 
// 	     pow(moments[1]/moments[0],x) * pow(moments[2]/moments[0],y)) * moments[0] * moments[0];

//   double s11 = sx * sx;
//   double s12 = rho * sx * sy;
//   double s21 = rho * sx * sy;
//   double s22 = sy * sy;

//   //Compute eigenvalues and eigenvectors of covariance matrix

//   double l1 = 0.5 * (s11 + s22 + sqrt((s11+s22)*(s11+s22) - 4*(s11*s22-s12*s21)));
//   double l2 = 0.5 * (s11 + s22 - sqrt((s11+s22)*(s11+s22) - 4*(s11*s22-s12*s21)));

//   if (l1 < 0)
//     l1 = 0.0;
//   if (l2 < 0)
//     l2 = 0.0;

//   double e11 = 1 / sqrt(1+((l1-s11)/s12)*((l1-s11)/s12));
//   double e12 = ((l1-s11)/s12) / sqrt(1+((l1-s11)/s12)*((l1-s11)/s12));
//   double e21 = 1 / sqrt(1+((l2-s11)/s12)*((l2-s11)/s12));
//   double e22 = ((l2-s11)/s12) / sqrt(1+((l2-s11)/s12)*((l2-s11)/s12));

//   //Gauss-Hermite Quadrature
//   int nodes = 4;

//   //Weights
//   double * w = new double[nodes];
//   if (nodes == 2)
//   {
//     w[0] = 8.862269254528e-1;
//     w[1] = 8.862269254528e-1;
//   }
//   else if (nodes == 4)
//   {
//     w[0] = 8.131283544725e-2;
//     w[1] = 8.049140900055e-1;
//     w[2] = 8.049140900055e-1;
//     w[3] = 8.131283544725e-2;
//   }

//   //Abscissas
//   double * a = new double[nodes];
//   if (nodes == 2)
//   {
//     a[0] = -0.707106781186548;
//     a[1] =  0.707106781186548;
//   }
//   else if (nodes == 4)
//   {
//     a[0] = -1.650680123885785;
//     a[1] = -0.524647623275290;
//     a[2] =  0.524647623275290;
//     a[3] =  1.650680123885785;
//   }

//   double source = 0;

//   for (int i = 0; i < nodes; i++)
//     for (int j = 0; j < nodes; j++)
//       for (int k = 0; k < nodes; k++)
// 	for (int l = 0; l < nodes; l++)
// 	  source += (1/(PI*PI)) * w[i] * w[j] * w[k] * w[l] * CoagFunc(sqrt(2*l1)*e11*a[i]+sqrt(2*l2)*e21*a[j]+ux,sqrt(2*l1)*e12*a[i]+sqrt(2*l2)*e22*a[j]+uy,sqrt(2*l1)*e11*a[k]+sqrt(2*l2)*e21*a[l]+ux,sqrt(2*l1)*e12*a[l]+sqrt(2*l2)*e22*a[k]+uy,x,y) * moments[0] * moments[0] * C;

//   delete [] w;
//   delete [] a;

//   return source;
// }

// double TSoot::CoagQuadLarge(double x, double y, double * moments, double temp)
// {
//   double C = GetC(temp);

//   double M00 = moments[0] - moments[6];
//   double M10 = moments[1] - 2.0*dimer_nbrC2*moments[6];
//   double M01 = moments[2] - pow(2.0*dimer_nbrC2, 2.0/3.0)*moments[6];

//   double M20 = moments[3] - pow(2.0*dimer_nbrC2, 2.0)*moments[6]; //M20
//   double M11 = moments[4] - pow(2.0*dimer_nbrC2, 5.0/3.0)*moments[6]; //M11
//   double M02 = moments[5] - pow(2.0*dimer_nbrC2, 4.0/3.0)*moments[6]; //M02

//   //double M20 = moments[3] - pow(2.0*dimer_nbrC2, 1.0/2.0)*moments[6]; //Mh0
//   //double M11 = moments[4] - pow(2.0*dimer_nbrC2, 5.0/6.0)*moments[6]; //Mhh
//   //double M02 = moments[5] - pow(2.0*dimer_nbrC2, 1.0/3.0)*moments[6]; //M0h

//   //Compute elements of covariance matrix
//   double ux = log(M10*M10/pow(M00*M00*M00*M20,1.0/2.0));
//   double uy = log(M01*M01/pow(M00*M00*M00*M02,1.0/2.0));
//   double sx = sqrt(MAX(log(M20*M00/(M10*M10)),0.0));
//   double sy = sqrt(MAX(log(M02*M00/(M01*M01)),0.0));
//   double rho = log(M11*M00/(M10*M01)) / (sx*sy);

//   //Check for delta function distribution
//   if (sx == 0 || sy == 0)
//   {
//     if (M00 < 1.0e-30)
//       return 0.0;
//     else
//     {
//       return C * pow(2.0,2.5) * pow(M10/M00,1.0/6.0) * ((1.0/2.0) * pow(2.0*M10/M00,x) * pow(2.0*M01/M00,y) - 
// 	     pow(M10/M00,x) * pow(M01/M00,y)) * M00 * M00;
//     }
//   }

//   double s11 = sx * sx;
//   double s12 = rho * sx * sy;
//   double s21 = rho * sx * sy;
//   double s22 = sy * sy;

//   //Compute eigenvalues and eigenvectors of covariance matrix
//   double l1 = 0.5 * (s11 + s22 + sqrt((s11+s22)*(s11+s22) - 4*(s11*s22-s12*s21)));
//   double l2 = 0.5 * (s11 + s22 - sqrt((s11+s22)*(s11+s22) - 4*(s11*s22-s12*s21)));

//   if (l1 < 0)
//     l1 = 0.0;
//   if (l2 < 0)
//     l2 = 0.0;

//   double e11 = 1 / sqrt(1+((l1-s11)/s12)*((l1-s11)/s12));
//   double e12 = ((l1-s11)/s12) / sqrt(1+((l1-s11)/s12)*((l1-s11)/s12));
//   double e21 = 1 / sqrt(1+((l2-s11)/s12)*((l2-s11)/s12));
//   double e22 = ((l2-s11)/s12) / sqrt(1+((l2-s11)/s12)*((l2-s11)/s12));

//   //Gauss-Hermite Quadrature
//   int nodes = 32;

//   //Weights
//   double * w = new double[nodes];
//   if (nodes == 1)
//   {
//     w[0] = 1.77245385091;
//   }
//   else if (nodes == 2)
//   {
//     w[0] = 8.862269254528e-1;
//     w[1] = 8.862269254528e-1;
//   }
//   else if (nodes == 4)
//   {
//     w[0] = 8.131283544725e-2;
//     w[1] = 8.049140900055e-1;
//     w[2] = 8.049140900055e-1;
//     w[3] = 8.131283544725e-2;
//   }
//   else if (nodes == 5)
//   {
//     w[0] = 1.995324205905e-2;
//     w[1] = 3.936193231522e-1;
//     w[2] = 9.453087204829e-1;
//     w[3] = 3.936193231522e-1;
//     w[4] = 1.995324205905e-2;
//   }
//   else if (nodes == 8)
//   {
//     w[0] = 1.996040722114e-4;
//     w[1] = 1.707798300741e-2;
//     w[2] = 2.078023258149e-1;
//     w[3] = 6.611470125582e-1;
//     w[4] = 6.611470125582e-1;
//     w[5] = 2.078023258149e-1;
//     w[6] = 1.707798300741e-2;
//     w[7] = 1.996040722114e-4;
//   }
//   else if (nodes == 16)
//   {
//     w[0] = 2.654807474011e-10;
//     w[1] = 2.320980844865e-7;
//     w[2] = 2.711860092538e-5;
//     w[3] = 9.322840086242e-4;
//     w[4] = 1.288031153551e-2;
//     w[5] = 8.381004139899e-2;
//     w[6] = 2.806474585285e-1;
//     w[7] = 5.079294790166e-1;
//     w[8] = w[7];
//     w[9] = w[6];
//     w[10] = w[5];
//     w[11] = w[4];
//     w[12] = w[3];
//     w[13] = w[2];
//     w[14] = w[1];
//     w[15] = w[0];
//   }
//   else if (nodes == 32)
//   {
//     w[15] = 3.75238352592802392864e-01;
//     w[14] = 2.77458142302529898131e-01;
//     w[13] = 1.51269734076642482578e-01;
//     w[12] = 6.04581309559126141860e-02;
//     w[11] = 1.75534288315734303030e-02;
//     w[10] = 3.65489032665442807915e-03;
//     w[ 9] = 5.36268365527972045989e-04;
//     w[ 8] = 5.41658406181998255789e-05;
//     w[ 7] = 3.65058512956237605727e-06;
//     w[ 6] = 1.57416779254559402923e-07;
//     w[ 5] = 4.09883216477089661816e-09;
//     w[ 4] = 5.93329146339663861478e-11;
//     w[ 3] = 4.21501021132644757306e-13;
//     w[ 2] = 1.19734401709284866582e-15;
//     w[ 1] = 9.23173653651829223381e-19;
//     w[ 0] = 7.31067642738416239302e-23;
//     w[16] = w[15];
//     w[17] = w[14];
//     w[18] = w[13];
//     w[19] = w[12];
//     w[20] = w[11];
//     w[21] = w[10];
//     w[22] = w[9];
//     w[23] = w[8];
//     w[24] = w[7];
//     w[25] = w[6];
//     w[26] = w[5];
//     w[27] = w[4];
//     w[28] = w[3];
//     w[29] = w[2];
//     w[30] = w[1];
//     w[31] = w[0];
//   }

//   //Abscissas
//   double * a = new double[nodes];
//   if (nodes == 1)
//   {
//     a[0] = 0.0;
//   }
//   else if (nodes == 2)
//   {
//     a[0] = -0.707106781186548;
//     a[1] =  0.707106781186548;
//   }
//   else if (nodes == 4)
//   {
//     a[0] = -1.650680123885785;
//     a[1] = -0.524647623275290;
//     a[2] =  0.524647623275290;
//     a[3] =  1.650680123885785;
//   }
//   else if (nodes == 5)
//   {
//     a[0] = -2.020182870456086;
//     a[1] = -0.958572464613819;
//     a[2] =  0.000000000000000;
//     a[3] =  0.958572464613819;
//     a[4] =  2.020182870456086;
//   }
//   else if (nodes == 8)
//   {
//     a[0] = -2.930637420257244;
//     a[1] = -1.981656756695843;
//     a[2] = -1.157193712446780;
//     a[3] = -0.381186990207322;
//     a[4] =  0.381186990207322;
//     a[5] =  1.157193712446780;
//     a[6] =  1.981656756695843;
//     a[7] =  2.930637420257244;
//   }
//   else if (nodes == 16)
//   {
//     a[0] = -4.68873893930582;
//     a[1] = -3.86944790486012;
//     a[2] = -3.17699916197996;
//     a[3] = -2.54620215784748;
//     a[4] = -1.95178799091625;
//     a[5] = -1.38025853919888;
//     a[6] = -0.82295144914466;
//     a[7] = -0.27348104613815;
//     a[8] = -a[7];
//     a[9] = -a[6];
//     a[10] = -a[5];
//     a[11] = -a[4];
//     a[12] = -a[3];
//     a[13] = -a[2];
//     a[14] = -a[1];
//     a[15] = -a[0];
//   }
//   else if (nodes == 32)
//   {
//     a[0] = -7.12581390983072757292e+00;
//     a[1] = -6.40949814926966041214e+00;
//     a[2] = -5.81222594951591383294e+00;
//     a[3] = -5.27555098651588012760e+00;
//     a[4] = -4.77716450350259639289e+00;
//     a[5] = -4.30554795335119844506e+00;
//     a[6] = -3.85375548547144464390e+00;
//     a[7] = -3.41716749281857073593e+00;
//     a[8] = -2.99249082500237420621e+00;
//     a[9] = -2.57724953773231745414e+00;
//     a[10] = -2.16949918360611217335e+00;
//     a[11] = -1.76765410946320160465e+00;
//     a[12] = -1.37037641095287183817e+00;
//     a[13] = -9.76500463589682838499e-01;
//     a[14] = -5.84978765435932448449e-01;
//     a[15] = -1.94840741569399326713e-01;
//     a[16] = -a[15];
//     a[17] = -a[14];
//     a[18] = -a[13];
//     a[19] = -a[12];
//     a[20] = -a[11];
//     a[21] = -a[10];
//     a[22] = -a[9];
//     a[23] = -a[8];
//     a[24] = -a[7];
//     a[25] = -a[6];
//     a[26] = -a[5];
//     a[27] = -a[4];
//     a[28] = -a[3];
//     a[29] = -a[2];
//     a[30] = -a[1];
//     a[31] = -a[0];
//   }

//   double source = 0;

//   for (int i = 0; i < nodes; i++)
//     for (int j = 0; j < nodes; j++)
//       for (int k = 0; k < nodes; k++)
// 	for (int l = 0; l < nodes; l++)
// 	  source += (1/(PI*PI)) * w[i] * w[j] * w[k] * w[l] * CoagFunc(sqrt(2*l1)*e11*a[i]+sqrt(2*l2)*e21*a[j]+ux,sqrt(2*l1)*e12*a[i]+sqrt(2*l2)*e22*a[j]+uy,sqrt(2*l1)*e11*a[k]+sqrt(2*l2)*e21*a[l]+ux,sqrt(2*l1)*e12*a[k]+sqrt(2*l2)*e22*a[l]+uy,x,y);

//   source *= M00 * M00 * C;
//   //source = sqrt(source)*M00 + pow(2.0*dimer_nbrC2,5.0/3.0)*moments[6];

//   delete [] w;
//   delete [] a;

//   return source;
// }

// double TSoot::CoagFunc(double x1, double y1, double x2, double y2, double i, double j)
// {
//   double eval;

//   double Df = 1.8;
//   double av = 1.0 - (2.0 / Df);
//   double as = (3.0 / Df) - 1.0;

//   //eval = exp(x1)*exp(y1)*exp(x2)*exp(y2);

//   eval = 0.5 * sqrt(exp(-x1) + exp(-x2)) * (exp(av*x1)*exp(as*y1)+exp(av*x2)*exp(as*y2)) * (exp(av*x1)*exp(as*y1)+exp(av*x2)*exp(as*y2)) * 
//     (pow(exp(x1)+exp(x2),i)*pow(exp(y1)+exp(y2),j)-exp(i*x1)*exp(j*y1)-exp(i*x2)*exp(j*y2));

//   return eval;
// }
