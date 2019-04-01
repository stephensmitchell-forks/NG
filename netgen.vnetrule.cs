namespace netgen
{

	public class vnetrule
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public vnetrule()
		{
		  name = new char[1];
		  name[0] = (char)0;
		  quality = 0;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void Dispose()
		{
		  // if (strlen(name)) 
		  name = null;
		  for (int i = 1; i <= freefaces.Size(); i++)
		  {
			freefaces.Elem(i) = null;
		  }
		  for (int i = 1; i <= freesets.Size(); i++)
		  {
			freesets.Elem(i) = null;
		  }
		  for (int i = 1; i <= freeedges.Size(); i++)
		  {
			freeedges.Elem(i) = null;
		  }
		  for (int i = 1; i <= freefaceinequ.Size(); i++)
		  {
			freefaceinequ.Elem(i) = null;
		  }
		  oldutofreezone = null;
		  oldutofreezonelimit = null;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int TestFlag(char flag)
		{
		  for (int i = 1; i <= flags.Size(); i++)
		  {
			if (flags.Get(i) == flag)
			{
				return 1;
			}
		  }
		  return 0;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void SetFreeZoneTransformation(Vector allp, int tolclass)
		{
		  int i;
		  int j;
		  // double nx, ny, nz, v1x, v1y, v1z, v2x, v2y, v2z;
		  double nl;
		  threeint ti;
		  int fs;
        
		  double lam1 = 1.0 / (2 * tolclass - 1);
		  double lam2 = 1 - lam1;
        
		  transfreezone.SetSize(freezone.Size());
        
		  int np = points.Size();
		  int nfp = freezone.Size();
		  Vector vp = new Vector(np);
		  Vector vfp1 = new Vector(nfp);
		  Vector vfp2 = new Vector(nfp);
        
        
		  for (i = 1; i <= 3; i++)
		  {
			  for (j = 1; j <= np; j++)
			  {
			vp(j - 1) = allp(i + 3 * j - 3 - 1);
			  }
        
			  oldutofreezone.Mult(vp, vfp1);
			  oldutofreezonelimit.Mult(vp, vfp2);
        
			  vfp1 *= lam1;
			  vfp1.Add.functorMethod(lam2, vfp2);
        
			  for (j = 1; j <= nfp; j++)
			  {
			transfreezone.Elem(j).X(i) = vfp1(j - 1);
			  }
		  }
        
		  // MARK(setfz2);
        
        
		  fzbox.SetPoint(transfreezone.Elem(1));
		  for (i = 2; i <= freezone.Size(); i++)
		  {
			fzbox.AddPoint(transfreezone.Elem(i));
		  }
        
        
		  // MARK(setfz3);
        
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  Array<threeint> freesetfaces = *freefaces.Get(fs);
			  DenseMatrix freesetinequ = *freefaceinequ.Get(fs);
        
			  for (i = 1; i <= freesetfaces.Size(); i++)
			  {
			  ti = freesetfaces.Get(i);
			  Point3d p1 = transfreezone.Get(ti.i1);
			  Point3d p2 = transfreezone.Get(ti.i2);
			  Point3d p3 = transfreezone.Get(ti.i3);
        
			  Vec3d v1 = new Vec3d(p1, p2);
			  Vec3d v2 = new Vec3d(p1, p3);
			  Vec3d n = new Vec3d();
			  Cross(v1, v2, n);
        
			  nl = n.Length();
        
			  if (nl < 1e-10)
			  {
				  freesetinequ.Set(1, 1, 0);
				  freesetinequ.Set(1, 2, 0);
				  freesetinequ.Set(1, 3, 0);
				  freesetinequ.Set(1, 4, -1);
			  }
			  else
			  {
				  //	      n /= nl;
        
				  freesetinequ.Set(i, 1, n.X() / nl);
				  freesetinequ.Set(i, 2, n.Y() / nl);
				  freesetinequ.Set(i, 3, n.Z() / nl);
				  freesetinequ.Set(i, 4, -(p1.X() * n.X() + p1.Y() * n.Y() + p1.Z() * n.Z()) / nl);
			  }
			  }
		  }
        
		  /*
		  (*testout) << "Transformed freezone: " << endl;
		  for (i = 1; i <= transfreezone.Size(); i++)
		    (*testout) << transfreezone.Get(i) << " ";
		  (*testout) << endl;
		  */
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int ConvexFreeZone()
		{
		  int i;
		  int j;
		  int k;
		  int fs;
        
		  // (*mycout) << "Convex free zone...\n";
        
		  int ret1 = 1;
		  // int ret2=1;
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  DenseMatrix freesetinequ = *freefaceinequ.Get(fs);
        
			  // const Array<int> & freeset = *freesets.Get(fs);
			  Array<twoint> freesetedges = *freeedges.Get(fs);
			  // const Array<threeint> & freesetfaces = *freefaces.Get(fs);
        
			  for (i = 1; i <= freesetedges.Size(); i++)
			  {
			  j = freesetedges.Get(i).i1; //triangle j with opposite point k
			  k = freesetedges.Get(i).i2;
        
			  if (freesetinequ.Get(j, 1) * transfreezone.Get(k).X() + freesetinequ.Get(j, 2) * transfreezone.Get(k).Y() + freesetinequ.Get(j, 3) * transfreezone.Get(k).Z() + freesetinequ.Get(j, 4) > 0)
			  {
				  ret1 = 0;
			  }
			  }
        
		  }
        
		  return ret1;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsInFreeZone(Point3d p)
		{
		  int i;
		  int fs;
		  char inthis;
        
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  inthis = 1;
			  Array<threeint> freesetfaces = *freefaces.Get(fs);
			  DenseMatrix freesetinequ = *freefaceinequ.Get(fs);
        
			  for (i = 1; i <= freesetfaces.Size() && inthis; i++)
			  {
			  if (freesetinequ.Get(i, 1) * p.X() + freesetinequ.Get(i, 2) * p.Y() + freesetinequ.Get(i, 3) * p.Z() + freesetinequ.Get(i, 4) > 0)
			  {
				inthis = 0;
			  }
			  }
        
			  if (inthis)
			  {
				  return 1;
			  }
		  }
        
		  return 0;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsTriangleInFreeZone(Point3d p1, Point3d p2, Point3d p3, Array<int> pi, int newone)
		{
		  int fs;
		  int infreeset;
		  int cannot = 0;
        
        
		  ArrayMem<int,3> pfi = new ArrayMem<int,3>(3);
		  ArrayMem<int,3> pfi2 = new ArrayMem<int,3>(3);
        
		  // convert from local index to freeset index
		  int i;
		  int j;
		  for (i = 1; i <= 3; i++)
		  {
			  pfi.Elem(i) = 0;
			  if (pi.Get(i))
			  {
			  for (j = 1; j <= freezonepi.Size(); j++)
			  {
				if (freezonepi.Get(j) == pi.Get(i))
				{
				  pfi.Elem(i) = j;
				}
			  }
			  }
		  }
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  Array<int> freeseti = *freesets.Get(fs);
			  for (i = 1; i <= 3; i++)
			  {
			  pfi2.Elem(i) = 0;
			  for (j = 1; j <= freeseti.Size(); j++)
			  {
				if (pfi.Get(i) == freeseti.Get(j))
				{
				  pfi2.Elem(i) = pfi.Get(i);
				}
			  }
			  }
        
			  infreeset = IsTriangleInFreeSet(p1, p2, p3, fs, pfi2, newone);
			  if (infreeset == 1)
			  {
				  return 1;
			  }
			  if (infreeset == -1)
			  {
				  cannot = -1;
			  }
		  }
        
		  return cannot;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsTriangleInFreeSet(Point3d p1, Point3d p2, Point3d p3, int fs, Array<int> pi, int newone)
		{
		  int i;
		  int ii;
		  Vec3d n = new Vec3d();
		  int allleft;
		  int allright;
		  int hos1;
		  int hos2;
		  int hos3;
		  int os1;
		  int os2;
		  int os3;
		  double hf;
		  double lam1;
		  double lam2;
		  double f;
		  double c1;
		  double c2;
		  double alpha;
		  double v1n;
		  double v2n;
		  double h11;
		  double h12;
		  double h22;
		  double dflam1;
		  double dflam2;
		  double lam1old;
		  double lam2old;
		  double fold;
		  double hpx;
		  double hpy;
		  double hpz;
		  double v1x;
		  double v1y;
		  double v1z;
		  double v2x;
		  double v2y;
		  double v2z;
		  int act1;
		  int act2;
		  int act3;
		  int it;
		  int cntout;
		  Array<int> activefaces = new Array<int>();
		  int isin;
        
        
		  // MARK(triinfz);
        
		  Array<threeint> freesetfaces = *freefaces.Get(fs);
		  DenseMatrix freesetinequ = *freefaceinequ.Get(fs);
        
        
		  int cnt = 0;
		  for (i = 1; i <= 3; i++)
		  {
			if (pi.Get(i))
			{
				cnt++;
			}
		  }
        
		  /*
		  (*testout) << "trig in free set : " << p1 << " - " << p2 << " - " << p3 << endl;
		  (*testout) << "common points: " << cnt << endl;
		  */
		  if (newone == 0)
		  {
			cnt = 0;
		  }
        
		  if (cnt == 1)
		  {
			  // MARK(triinfz1);
        
			  int upi = 0;
			  int lpiu = 0;
			  for (i = 1; i <= 3; i++)
			  {
			if (pi.Get(i))
			{
				upi = i;
				lpiu = pi.Get(i);
			}
			  }
        
			  Vec3d v1 = new Vec3d();
			  Vec3d v2 = new Vec3d();
			  switch (upi)
			  {
			case 1:
			{
				v1 = p2 - p1;
				v2 = p3 - p1;
				break;
			}
			case 2:
			{
				v1 = p3 - p2;
				v2 = p1 - p2;
				break;
			}
			case 3:
			{
				v1 = p1 - p3;
				v2 = p2 - p3;
				break;
			}
			  }
        
			  v1 /= v1.Length();
			  v2 /= v2.Length();
			  Cross(v1, v2, n);
			  n /= n.Length();
        
			  //      (*testout) << "Test new: " << endl;
			  for (i = 1; i <= freesetfaces.Size(); i++)
			  {
			  if ((freesetfaces.Get(i).i1 == lpiu) || (freesetfaces.Get(i).i2 == lpiu) || (freesetfaces.Get(i).i3 == lpiu))
			  {
				  // freeface has point
        
        
				  Vec3d a = new Vec3d(freesetinequ.Get(i, 1), freesetinequ.Get(i, 2), freesetinequ.Get(i, 3));
        
				  //	      if (1 - fabs (a * n) < 1e-8 ) 
				  //		continue;
        
				  Vec3d an = new Vec3d();
				  Cross(a, n, an);
				  double lan = an.Length();
				  if (lan < 1e-10)
				  {
				continue;
				  }
        
				  an /= lan;
        
				  int out1 = (a * v1) > 0;
				  int out2 = (a * v2) > 0;
				  //	      (*testout) << "out1, out2 = " << out1 << ", " << out2 << endl;
				  if (out1 != 0 && out2 != 0)
				  {
				return 0;
				  }
        
				  if (out1 == 0 && out2 == 0)
				  {
				continue;
				  }
        
        
				  //	      if ( ( (an * v1) < 0) &&  ( (an * v2) < 0) )   // falsch !!!!
				  //		an *= -1;
        
				  // solve  an = lam1 v1 + lam2 v2
				  double vii11 = v1 * v1;
				  double vii12 = v1 * v2;
				  double vii22 = v2 * v2;
				  double det = vii11 * vii22 - vii12 * vii12;
				  if (ngsimd.GlobalMembers.fabs(det) < 1e-10)
				  {
				continue;
				  }
				  double rs1 = an * v1;
				  double rs2 = an * v2;
        
				  double lambda1 = rs1 * vii22 - rs2 * vii12;
				  double lambda2 = rs2 * vii11 - rs1 * vii12;
        
				  if (ngsimd.GlobalMembers.fabs(lambda1) > ngsimd.GlobalMembers.fabs(lambda2))
				  {
				  if (lambda1 < 0)
				  {
					an *= -1;
				  }
				  }
				  else
				  {
				  if (lambda2 < 0)
				  {
					an *= -1;
				  }
				  }
        
        
				  if (lambda1 * lambda2 < 0 && false)
				  {
				  if (ngsimd.GlobalMembers.fabs(lambda1) > 1e-14 && ngsimd.GlobalMembers.fabs(lambda2) > 1e-14)
				  {
					  //		      (*mycout) << "lambda1 lambda2 < 0" << endl;
					  (*testout) << "lambdai different" << "\n";
					  (*testout) << "v1 = " << v1 << "\n";
					  (*testout) << "v2 = " << v2 << "\n";
					  (*testout) << "n = " << n << "\n";
					  (*testout) << "a = " << a << "\n";
					  (*testout) << "an = " << an << "\n";
					  (*testout) << "a * v1 = " << (a * v1) << "\n";
					  (*testout) << "a * v2 = " << (a * v2) << "\n";
					  (*testout) << "an * v1 = " << (an * v1) << "\n";
					  (*testout) << "an * v2 = " << (an * v2) << "\n";
        
					  (*testout) << "vii = " << vii11 << ", " << vii12 << ", " << vii22 << "\n";
					  (*testout) << "lambdai = " << lambda1 << ", " << lambda2 << "\n";
					  (*testout) << "rs = " << rs1 << ", " << rs2 << "\n";
					  continue;
				  }
				  }
        
				  if (out1 != 0)
				  {
				v1.CopyFrom(an);
				  }
				  else
				  {
				v2.CopyFrom(an);
				  }
			  }
			  }
        
			  return 1;
        
			  /*
			  (*testout) << "overlap trig " << p1 << p2 << p3 << endl;
			  (*testout) << "upi = " << upi << endl;
			  (*testout) << "v1 = " << v1 << " v2 = " << v2 << endl;
			  */
        
			  switch (upi)
			  {
			case 1:
			{
				v1 = p2 - p1;
				v2 = p3 - p1;
				break;
			}
			case 2:
			{
				v1 = p3 - p2;
				v2 = p1 - p2;
				break;
			}
			case 3:
			{
				v1 = p1 - p3;
				v2 = p2 - p3;
				break;
			}
			  }
        
			  v1 /= v1.Length();
			  v2 /= v2.Length();
			  Cross(v1, v2, n);
			  n /= n.Length();
        
			  //      (*testout) << "orig v1, v2 = " << v1 << ", " << v2 << endl;
        
        
			  for (i = 1; i <= freesetfaces.Size(); i++)
			  {
			  if ((freesetfaces.Get(i).i1 == lpiu) || (freesetfaces.Get(i).i2 == lpiu) || (freesetfaces.Get(i).i3 == lpiu))
			  {
				  /*
				  (*testout) << "v1, v2, now = " << v1 << ", " << v2 << endl;
		
				  // freeface has point
				  (*testout) << "freesetface: "
					 << freesetfaces.Get(i).i1 << " "
					 << freesetfaces.Get(i).i2 << " "
					 << freesetfaces.Get(i).i3 << " ";
				  */
        
				  Vec3d a = new Vec3d(freesetinequ.Get(i, 1), freesetinequ.Get(i, 2), freesetinequ.Get(i, 3));
				  //	      (*testout) << "a = " <<  a << endl;
        
        
				  Vec3d an = new Vec3d();
				  Cross(a, n, an);
				  double lan = an.Length();
        
				  //	      (*testout) << "an = " << an << endl;
        
				  if (lan < 1e-10)
				  {
				continue;
				  }
        
				  an /= lan;
        
				  //	      (*testout) << "a*v1 = " << (a*v1) << " a*v2 = " << (a*v2) << endl;
        
				  int out1 = (a * v1) > 0;
				  // int out2 = (a * v2) > 0;
        
        
				  //	      (*testout) << "out1, 2 = " << out1 << ", " << out2 << endl;
        
        
				  double vii11 = v1 * v1;
				  double vii12 = v1 * v2;
				  double vii22 = v2 * v2;
				  double det = vii11 * vii22 - vii12 * vii12;
				  if (ngsimd.GlobalMembers.fabs(det) < 1e-10)
				  {
				continue;
				  }
				  double rs1 = an * v1;
				  double rs2 = an * v2;
        
				  double lambda1 = rs1 * vii22 - rs2 * vii12;
				  double lambda2 = rs2 * vii11 - rs1 * vii12;
        
				  //	      (*testout) << "lambda1, lambda2 = " << lambda1 << ", " << lambda2 << endl;
        
        
				  if (ngsimd.GlobalMembers.fabs(lambda1) > ngsimd.GlobalMembers.fabs(lambda2))
				  {
				  if (lambda1 < 0)
				  {
					an *= -1;
				  }
				  }
				  else
				  {
				  if (lambda2 < 0)
				  {
					an *= -1;
				  }
				  }
        
        
				  if (lambda1 * lambda2 < 0)
				  {
				  if (ngsimd.GlobalMembers.fabs(lambda1) > 1e-14 && ngsimd.GlobalMembers.fabs(lambda2) > 1e-14)
				  {
					  //		      (*mycout) << "lambda1 lambda2 < 0" << endl;
					  (*testout) << "lambdai different" << "\n";
					  (*testout) << "v1 = " << v1 << "\n";
					  (*testout) << "v2 = " << v2 << "\n";
					  (*testout) << "n = " << n << "\n";
					  (*testout) << "a = " << a << "\n";
					  (*testout) << "an = " << an << "\n";
					  (*testout) << "a * v1 = " << (a * v1) << "\n";
					  (*testout) << "a * v2 = " << (a * v2) << "\n";
					  (*testout) << "an * v1 = " << (an * v1) << "\n";
					  (*testout) << "an * v2 = " << (an * v2) << "\n";
        
					  (*testout) << "vii = " << vii11 << ", " << vii12 << ", " << vii22 << "\n";
					  (*testout) << "lambdai = " << lambda1 << ", " << lambda2 << "\n";
					  (*testout) << "rs = " << rs1 << ", " << rs2 << "\n";
					  continue;
				  }
				  }
        
				  if (out1 != 0)
				  {
				v1.CopyFrom(an);
				  }
				  else
				  {
				v2.CopyFrom(an);
				  }
        
        
        
			  }
			  }
        
			  return 1;
		  }
        
        
        
		  if (cnt == 2)
		  {
			  //      (*testout) << "tripoitns: " << p1 << " " << p2 << " " << p3 << endl;
        
			  // MARK(triinfz2);
        
			  int pi1 = 0;
			  int pi2 = 0;
			  int pi3 = 0;
			  Vec3d a1 = new Vec3d(); // outer normals
			  Vec3d a2 = new Vec3d();
			  Vec3d trivec = new Vec3d(); // vector from common edge to third point of triangle
			  for (i = 1; i <= 3; i++)
			  {
			if (pi.Get(i))
			{
				pi2 = pi1;
				pi1 = pi.Get(i);
			}
			else
			{
			  pi3 = i;
			}
			  }
        
			  switch (pi3)
			  {
			case 1:
				trivec = (p1 - p2);
				break;
			case 2:
				trivec = (p2 - p3);
				break;
			case 3:
				trivec = (p3 - p2);
				break;
			  }
        
			  Array<int> lpi = new Array<int>(freezonepi.Size());
			  for (i = 1; i <= lpi.Size(); i++)
			  {
			lpi.Elem(i) = 0;
			  }
			  lpi.Elem(pi1) = 1;
			  lpi.Elem(pi2) = 1;
        
			  int ff1 = 0;
			  int ff2 = 0;
			  for (i = 1; i <= freesetfaces.Size(); i++)
			  {
			  if (lpi.Get(freesetfaces.Get(i).i1) + lpi.Get(freesetfaces.Get(i).i2) + lpi.Get(freesetfaces.Get(i).i3) == 2)
			  {
				  ff2 = ff1;
				  ff1 = i;
			  }
			  }
        
			  if (ff2 == 0)
			  {
			return 1;
			  }
        
			  a1 = new Vec3d(freesetinequ.Get(ff1, 1), freesetinequ.Get(ff1, 2), freesetinequ.Get(ff1, 3));
			  a2 = new Vec3d(freesetinequ.Get(ff2, 1), freesetinequ.Get(ff2, 2), freesetinequ.Get(ff2, 3));
        
			  if (((a1 * trivec) > 0) || ((a2 * trivec) > 0))
			  {
			return 0;
			  }
        
			  return 1;
		  }
        
        
		  if (cnt == 3)
		  {
			  // MARK(triinfz3);  
        
			  Array<int> lpi = new Array<int>(freezonepi.Size());
			  for (i = 1; i <= lpi.Size(); i++)
			  {
			lpi.Elem(i) = 0;
			  }
        
			  for (i = 1; i <= 3; i++)
			  {
			lpi.Elem(pi.Get(i)) = 1;
			  }
        
			  for (i = 1; i <= freesetfaces.Size(); i++)
			  {
			  if (lpi.Get(freesetfaces.Get(i).i1) + lpi.Get(freesetfaces.Get(i).i2) + lpi.Get(freesetfaces.Get(i).i3) == 3)
			  {
				  return 0;
			  }
			  }
			  return 1;
		  }
        
		  // MARK(triinfz0);  
        
        
		  os1 = os2 = os3 = 0;
		  activefaces.SetSize(0);
        
		  // is point inside ?
        
		  for (i = 1; i <= freesetfaces.Size(); i++)
		  {
			  hos1 = (int)freesetinequ.Get(i, 1) * p1.X() + freesetinequ.Get(i, 2) * p1.Y() + freesetinequ.Get(i, 3) * p1.Z() + freesetinequ.Get(i, 4) > -1E-5;
        
			  hos2 = (int)freesetinequ.Get(i, 1) * p2.X() + freesetinequ.Get(i, 2) * p2.Y() + freesetinequ.Get(i, 3) * p2.Z() + freesetinequ.Get(i, 4) > -1E-5;
        
			  hos3 = (int)freesetinequ.Get(i, 1) * p3.X() + freesetinequ.Get(i, 2) * p3.Y() + freesetinequ.Get(i, 3) * p3.Z() + freesetinequ.Get(i, 4) > -1E-5;
        
			  if (hos1 != 0 && hos2 != 0 && hos3 != 0)
			  {
				  return 0;
			  }
        
			  if (hos1 != 0)
			  {
				  os1 = 1;
			  }
			  if (hos2 != 0)
			  {
				  os2 = 1;
			  }
			  if (hos3 != 0)
			  {
				  os3 = 1;
			  }
        
			  if (hos1 != 0 || hos2 != 0 || hos3 != 0)
			  {
				  activefaces.Append(i);
			  }
		  }
        
		  if (os1 == 0 || os2 == 0 || os3 == 0)
		  {
			  return 1;
		  }
        
		  v1x = p2.X() - p1.X();
		  v1y = p2.Y() - p1.Y();
		  v1z = p2.Z() - p1.Z();
        
		  v2x = p3.X() - p1.X();
		  v2y = p3.Y() - p1.Y();
		  v2z = p3.Z() - p1.Z();
        
		  n.X() = v1y * v2z - v1z * v2y;
		  n.Y() = v1z * v2x - v1x * v2z;
		  n.Z() = v1x * v2y - v1y * v2x;
		  n /= n.Length();
        
		  allleft = allright = 1;
		  for (i = 1; i <= transfreezone.Size() && (allleft || allright); i++)
		  {
			  Point3d p = transfreezone.Get(i);
			  float scal = (p.X() - p1.X()) * n.X() + (p.Y() - p1.Y()) * n.Y() + (p.Z() - p1.Z()) * n.Z();
        
			  if (scal > 1E-8F)
			  {
				  allleft = 0;
			  }
			  if (scal < -1E-8F)
			  {
				  allright = 0;
			  }
		  }
        
		  if (allleft != 0 || allright != 0)
		  {
			  return 0;
		  }
        
        
		  lam1old = lam2old = lam1 = lam2 = 1.0 / 3.0;
        
        
		  //  testout << endl << endl << "Start minimizing" << endl;
        
		  it = 0;
		  int minit;
		  minit = 1000;
		  fold = 1E10;
        
        
        
		  while (true)
		  {
			  it++;
        
			  if (it > 1000)
			  {
				  return -1;
			  }
        
			  if (lam1 < 0)
			  {
				  lam1 = 0;
			  }
			  if (lam2 < 0)
			  {
				  lam2 = 0;
			  }
			  if (lam1 + lam2 > 1)
			  {
				  lam1 = 1 - lam2;
			  }
        
			  if (it > minit)
			  {
			  (*testout) << "it = " << it << "\n";
			  (*testout) << "lam1/2 = " << lam1 << "  " << lam2 << "\n";
			  }
        
			  hpx = p1.X() + lam1 * v1x + lam2 * v2x;
			  hpy = p1.Y() + lam1 * v1y + lam2 * v2y;
			  hpz = p1.Z() + lam1 * v1z + lam2 * v2z;
        
			  f = 0;
        
			  h11 = h12 = h22 = dflam1 = dflam2 = 0;
			  cntout = 0;
        
			  isin = 1;
        
			  for (i = 1; i <= activefaces.Size(); i++)
			  {
			  ii = activefaces.Get(i);
        
			  hf = freesetinequ.Get(ii, 1) * hpx + freesetinequ.Get(ii, 2) * hpy + freesetinequ.Get(ii, 3) * hpz + freesetinequ.Get(ii, 4);
        
			  if (hf > -1E-7)
			  {
				  isin = 0;
			  }
        
			  hf += 1E-4;
			  if (hf > 0)
			  {
				  f += hf * hf;
        
				  v1n = freesetinequ.Get(ii, 1) * v1x + freesetinequ.Get(ii, 2) * v1y + freesetinequ.Get(ii, 3) * v1z;
				  v2n = freesetinequ.Get(ii, 1) * v2x + freesetinequ.Get(ii, 2) * v2y + freesetinequ.Get(ii, 3) * v2z;
        
				  h11 += 2 * v1n * v1n;
				  h12 += 2 * v1n * v2n;
				  h22 += 2 * v2n * v2n;
				  dflam1 += 2 * hf * v1n;
				  dflam2 += 2 * hf * v2n;
				  cntout++;
			  }
			  }
        
			  if (isin != 0)
			  {
				  return 1;
			  }
        
			  if (it > minit)
			  {
			  (*testout) << "f = " << f << "  dfdlam = " << dflam1 << "  " << dflam2 << "\n";
			  (*testout) << "h = " << h11 << "  " << h12 << "  " << h22 << "\n";
			  (*testout) << "active: " << cntout << "\n";
			  (*testout) << "lam1-lam1old = " << (lam1 - lam1old) << "\n";
			  (*testout) << "lam2-lam2old = " << (lam2 - lam2old) << "\n";
			  }
        
        
			  if (f >= fold)
			  {
			  lam1 = 0.100000000000000 * lam1 + 0.9000000000000000 * lam1old;
			  lam2 = 0.100000000000000 * lam2 + 0.9000000000000000 * lam2old;
			  }
			  else
			  {
			  lam1old = lam1;
			  lam2old = lam2;
			  fold = f;
        
        
			  if (f < 1E-9)
			  {
				  return 1;
			  }
        
			  h11 += 1E-10;
			  h22 += 1E-10;
			  c1 = - (h22 * dflam1 - h12 * dflam2) / (h11 * h22 - h12 * h12);
			  c2 = - (-h12 * dflam1 + h11 * dflam2) / (h11 * h22 - h12 * h12);
			  alpha = 1;
        
        
			  if (it > minit)
			  {
				(*testout) << "c1/2 = " << c1 << "  " << c2 << "\n";
			  }
        
			  act1 = lam1 <= 1E-6 && c1 <= 0;
			  act2 = lam2 <= 1E-6 && c2 <= 0;
			  act3 = lam1 + lam2 >= 1 - 1E-6 && c1 + c2 >= 0;
        
			  if (it > minit)
			  {
				(*testout) << "act1,2,3 = " << act1 << act2 << act3 << "\n";
			  }
        
			  if ((act1 && act2) || (act1 && act3) || (act2 && act3))
			  {
				  return 0;
			  }
        
			  if (act1 != 0)
			  {
				  c1 = 0;
				  c2 = - dflam2 / h22;
			  }
        
			  if (act2 != 0)
			  {
				  c1 = - dflam1 / h11;
				  c2 = 0;
			  }
        
			  if (act3 != 0)
			  {
				  c1 = - (dflam1 - dflam2) / (h11 + h22 - 2 * h12);
				  c2 = -c1;
			  }
        
			  if (it > minit)
			  {
				(*testout) << "c1/2 now = " << c1 << "  " << c2 << "\n";
			  }
        
        
			  if (f > 100 * ngsimd.GlobalMembers.sqrt(sqr(c1) + sqr(c2)))
			  {
				  return 0;
			  }
        
        
			  if (lam1 + alpha * c1 < 0 && act1 == 0)
			  {
				alpha = -lam1 / c1;
			  }
			  if (lam2 + alpha * c2 < 0 && act2 == 0)
			  {
				alpha = -lam2 / c2;
			  }
			  if (lam1 + lam2 + alpha * (c1 + c2) > 1 && act3 == 0)
			  {
				alpha = (1 - lam1 - lam2) / (c1 + c2);
			  }
        
			  if (it > minit)
			  {
				(*testout) << "alpha = " << alpha << "\n";
			  }
        
			  lam1 += alpha * c1;
			  lam2 += alpha * c2;
			  }
		  }
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsQuadInFreeZone(Point3d p1, Point3d p2, Point3d p3, Point3d p4, Array<int> pi, int newone)
		{
		  int fs;
		  int infreeset;
		  int cannot = 0;
        
        
		  ArrayMem<int,4> pfi = new ArrayMem<int,4>(4);
		  ArrayMem<int,4> pfi2 = new ArrayMem<int,4>(4);
        
		  // convert from local index to freeset index
		  int i;
		  int j;
		  for (i = 1; i <= 4; i++)
		  {
			  pfi.Elem(i) = 0;
			  if (pi.Get(i))
			  {
			  for (j = 1; j <= freezonepi.Size(); j++)
			  {
				if (freezonepi.Get(j) == pi.Get(i))
				{
				  pfi.Elem(i) = j;
				}
			  }
			  }
		  }
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  Array<int> freeseti = *freesets.Get(fs);
			  for (i = 1; i <= 4; i++)
			  {
			  pfi2.Elem(i) = 0;
			  for (j = 1; j <= freeseti.Size(); j++)
			  {
				if (pfi.Get(i) == freeseti.Get(j))
				{
				  pfi2.Elem(i) = pfi.Get(i);
				}
			  }
			  }
        
			  infreeset = IsQuadInFreeSet(p1, p2, p3, p4, fs, pfi2, newone);
			  if (infreeset == 1)
			  {
				  return 1;
			  }
			  if (infreeset == -1)
			  {
				  cannot = -1;
			  }
		  }
        
		  return cannot;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsQuadInFreeSet(Point3d p1, Point3d p2, Point3d p3, Point3d p4, int fs, Array<int> pi, int newone)
		{
		  int i;
        
		  int cnt = 0;
		  for (i = 1; i <= 4; i++)
		  {
			if (pi.Get(i))
			{
				cnt++;
			}
		  }
        
		  /*
		  (*testout) << "test quad in freeset: " << p1 << " - " << p2 << " - " << p3 << " - " << p4 << endl;
		  (*testout) << "pi = ";
		  for (i = 1; i <= pi.Size(); i++)
		    (*testout) << pi.Get(i) << " ";
		  (*testout) << endl;
		  (*testout) << "cnt = " << cnt  << endl;
		  */
		  if (cnt == 4)
		  {
			  return 1;
		  }
        
		  if (cnt == 3)
		  {
			  return 1;
		  }
        
		  ArrayMem<int,3> pi3 = new ArrayMem<int,3>(3);
		  int res;
        
		  pi3.Elem(1) = pi.Get(1);
		  pi3.Elem(2) = pi.Get(2);
		  pi3.Elem(3) = pi.Get(3);
		  res = IsTriangleInFreeSet(p1, p2, p3, fs, pi3, newone);
		  if (res != 0)
		  {
			  return res;
		  }
        
        
		  pi3.Elem(1) = pi.Get(2);
		  pi3.Elem(2) = pi.Get(3);
		  pi3.Elem(3) = pi.Get(4);
		  res = IsTriangleInFreeSet(p2, p3, p4, fs, pi3, newone);
		  if (res != 0)
		  {
			  return res;
		  }
        
		  pi3.Elem(1) = pi.Get(3);
		  pi3.Elem(2) = pi.Get(4);
		  pi3.Elem(3) = pi.Get(1);
		  res = IsTriangleInFreeSet(p3, p4, p1, fs, pi3, newone);
		  if (res != 0)
		  {
			  return res;
		  }
        
		  pi3.Elem(1) = pi.Get(4);
		  pi3.Elem(2) = pi.Get(1);
		  pi3.Elem(3) = pi.Get(2);
		  res = IsTriangleInFreeSet(p4, p1, p2, fs, pi3, newone);
		  return res;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public float CalcPointDist(int pi, Point3d p)
		{
		  float dx = p.X() - points.Get(pi).X();
		  float dy = p.Y() - points.Get(pi).Y();
		  float dz = p.Z() - points.Get(pi).Z();
        
		  //  const threefloat * tf = &tolerances.Get(pi);
		  //  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
		  return tolerances.Get(pi) * (dx * dx + dy * dy + dz * dz);
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int TestOk()
		{
		  Array<int> cntpused = new Array<int>(points.Size());
		  Array<int> edge1 = new Array<int>();
		  Array<int> edge2 = new Array<int>();
		  Array<int> delf = new Array<int>(faces.Size());
		  int i;
		  int j;
		  int k;
		  int pi1;
		  int pi2;
		  int found;
        
		  for (i = 1; i <= cntpused.Size(); i++)
		  {
			cntpused.Elem(i) = 0;
		  }
		  for (i = 1; i <= faces.Size(); i++)
		  {
			delf.Elem(i) = 0;
		  }
		  for (i = 1; i <= delfaces.Size(); i++)
		  {
			delf.Elem(delfaces.Get(i)) = 1;
		  }
        
        
		  for (i = 1; i <= faces.Size(); i++)
		  {
			if (delf.Get(i) || i > noldf)
			{
			  for (j = 1; j <= faces.Get(i).GetNP(); j++)
			  {
				cntpused.Elem(faces.Get(i).PNum(j))++;
			  }
			}
		  }
        
		  for (i = 1; i <= cntpused.Size(); i++)
		  {
			if (cntpused.Get(i) > 0 && cntpused.Get(i) < 2)
			{
			return 0;
			}
		  }
        
        
		  //  (*testout) << endl;
		  for (i = 1; i <= faces.Size(); i++)
		  {
			  //      (*testout) << "face " << i << endl;
			  for (j = 1; j <= faces.Get(i).GetNP(); j++)
			  {
			  pi1 = 0;
			  pi2 = 0;
			  if (delf.Get(i))
			  {
				  pi1 = faces.Get(i).PNumMod(j);
				  pi2 = faces.Get(i).PNumMod(j + 1);
			  }
			  if (i > noldf)
			  {
				  pi1 = faces.Get(i).PNumMod(j + 1);
				  pi2 = faces.Get(i).PNumMod(j);
			  }
        
			  found = 0;
			  if (pi1 != 0)
			  {
				  for (k = 1; k <= edge1.Size(); k++)
				  {
				if (edge1.Get(k) == pi1 && edge2.Get(k) == pi2)
				{
					found = 1;
					edge1.DeleteElement(k);
					edge2.DeleteElement(k);
					k--;
					//		    (*testout) << "Del edge " << pi1 << "-" << pi2 << endl;
				}
				  }
				  if (found == 0)
				  {
				  edge1.Append(pi2);
				  edge2.Append(pi1);
				  //		  (*testout) << "Add edge " << pi1 << "-" << pi2 << endl;
				  }
			  }
			  }
		  }
        
        
		  if (edge1.Size() > 0)
		  {
			  return 0;
		  }
        
		  /*
		    cntpused.SetSize(freezone.Size());
		    for (i = 1; i <= cntpused.Size(); i++)
		    cntpused[i] = 0;
		
		    for (i = 1; i <= freefaces.Size(); i++)
		    {
		    cntpused[freefaces[i].i1]++;
		    cntpused[freefaces[i].i2]++;
		    cntpused[freefaces[i].i3]++;
		    }
		
		    for (i = 1; i <= cntpused.Size(); i++)
		    if (cntpused[i] < 3)
		    {
		    (*mycout) << "Fall 3" << endl;
		    return 0;
		    }
		
		
		
		    for (i = 1; i <= freefaces.Size(); i++)
		    {
		    for (j = 1; j <= 3; j++)
		    {
		    if (j == 1)
		    {
		    pi1 = freefaces[i].i1;
		    pi2 = freefaces[i].i2;
		    }
		    if (j == 2)
		    {
		    pi1 = freefaces[i].i2;
		    pi2 = freefaces[i].i3;
		    }
		    if (j == 3)
		    {
		    pi1 = freefaces[i].i3;
		    pi2 = freefaces[i].i1;
		    }
		
		    found = 0;
		    for (k = 1; k <= edge1.Size(); k++)
		    if (edge1[k] == pi1 && edge2[k] == pi2)
		    {
		    found = 1;
		    edge1.DeleteElement(k);
		    edge2.DeleteElement(k);
		    k--;
		    }
		
		    if (!found)
		    {
		    edge1.Append (pi2);
		    edge2.Append (pi1);
		    }
		    }
		    }
		
		    if (edge1.Size() > 0)
		    {
		    (*mycout) << "Fall 4" << endl;
		    return 0;
		    }
		    */
		  return 1;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsDelFace(int fn)
		{
		  int i;
		  for (i = 1; i <= GetNDelF(); i++)
		  {
			if (GetDelFace(i) == fn)
			{
				return 1;
			}
		  }
		  return 0;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int NeighbourTrianglePoint(threeint t1, threeint t2)
		{
		  Array<int> tr1 = new Array<int>(3);
		  Array<int> tr2 = new Array<int>(3);
		  tr1.Elem(1) = t1.i1;
		  tr1.Elem(2) = t1.i2;
		  tr1.Elem(3) = t1.i3;
		  tr2.Elem(1) = t2.i1;
		  tr2.Elem(2) = t2.i2;
		  tr2.Elem(3) = t2.i3;
        
        
		  int ret = 0;
        
		  for (int i = 1; i <= 3; i++)
		  {
			  for (int j = 1; j <= 3; j++)
			  {
			  if ((tr1.Get(i) == tr2.Get(j) && tr1.Get((i % 3) + 1) == tr2.Get((j % 3) + 1)) || (tr1.Get(i) == tr2.Get((j % 3) + 1) && tr1.Get((i % 3) + 1) == tr2.Get(j)))
			  {
					ret = tr2.Get((j + 1) % 3 + 1);
			  }
			  }
		  }
        
		  return ret;
        
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void LoadRule(istream ist)
		{
		  string buf = new string(new char[256]);
		  char ch;
		  char ok;
		  Point3d p = new Point3d();
		  Element2d face = new Element2d(ELEMENT_TYPE.TRIG);
		  int i;
		  int j;
		  int i1;
		  int i2;
		  int i3;
		  int fs;
		  int ii;
		  int ii1;
		  int ii2;
		  int ii3;
		  twoint edge = new twoint();
		  DenseMatrix tempoldutonewu = new DenseMatrix(30, 20);
		  DenseMatrix tempoldutofreezone = new DenseMatrix(30, 20);
		  DenseMatrix tempoldutofreezonelimit = new DenseMatrix(30, 20);
		  DenseMatrix tfz = new DenseMatrix(20, 20);
		  DenseMatrix tfzl = new DenseMatrix(20, 20);
        
		  tempoldutonewu = 0;
		  tempoldutofreezone = 0;
		  tfz = 0;
		  tfzl = 0;
        
        
		  noldp = 0;
		  noldf = 0;
        
		  ist.get(buf, sizeof(char), '"');
		  ist.get(ch);
		  ist.get(buf, sizeof(char), '"');
		  ist.get(ch);
        
		  name = null;
		  name = new char[buf.Length + 1];
		  name = buf;
		  //  (*mycout) << "Rule " << name << " found." << endl;
        
		  do
		  {
			  ist >> buf;
        
			  if (string.Compare(buf, "quality") == 0)
        
			  {
			  ist >> quality;
			  }
        
			  else if (string.Compare(buf, "flags") == 0)
			  {
			  ist >> ch;
			  while (ch != ';')
			  {
				  flags.Append(ch);
				  ist >> ch;
			  }
			  }
        
			  else if (string.Compare(buf, "mappoints") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ','
				  ist >> p.Z();
				  ist >> ch; // ')'
        
				  points.Append(p);
				  noldp++;
        
				  tolerances.SetSize(noldp);
				  tolerances.Elem(noldp) = 1;
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  ist >> tolerances.Elem(noldp);
					  ist >> ch; // '}'
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
        
			  else if (string.Compare(buf, "mapfaces") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  face.SetType(ELEMENT_TYPE.TRIG);
				  ist >> face.PNum(1);
				  ist >> ch; // ','
				  ist >> face.PNum(2);
				  ist >> ch; // ','
				  ist >> face.PNum(3);
				  ist >> ch; // ')' or ','
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  face.SetType(ELEMENT_TYPE.QUAD);
				  ist >> face.PNum(4);
				  ist >> ch; // ')'
				  }
				  faces.Append(face);
				  noldf++;
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == 'd')
				  {
					  delfaces.Append(noldf);
					  ist >> ch; // 'e'
					  ist >> ch; // 'l'
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "mapedges") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> edge.i1;
				  ist >> ch; // ','
				  ist >> edge.i2;
				  ist >> ch; // ')'
        
				  edges.Append(edge);
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
        
			  else if (string.Compare(buf, "newpoints") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ','
				  ist >> p.Z();
				  ist >> ch; // ')'
        
				  points.Append(p);
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  LoadVMatrixLine(ist, tempoldutonewu.functorMethod, 3 * (points.Size() - noldp) - 2);
        
					  ist >> ch; // '{'
					  LoadVMatrixLine(ist, tempoldutonewu.functorMethod, 3 * (points.Size() - noldp) - 1);
        
					  ist >> ch; // '{'
					  LoadVMatrixLine(ist, tempoldutonewu.functorMethod, 3 * (points.Size() - noldp));
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "newfaces") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  face.SetType(ELEMENT_TYPE.TRIG);
				  ist >> face.PNum(1);
				  ist >> ch; // ','
				  ist >> face.PNum(2);
				  ist >> ch; // ','
				  ist >> face.PNum(3);
				  ist >> ch; // ')' or ','
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  face.SetType(ELEMENT_TYPE.QUAD);
				  ist >> face.PNum(4);
				  ist >> ch; // ')'
				  }
				  faces.Append(face);
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "freezone") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ','
				  ist >> p.Z();
				  ist >> ch; // ')'
        
				  freezone.Append(p);
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  LoadVMatrixLine(ist, tempoldutofreezone.functorMethod, 3 * freezone.Size() - 2);
        
					  ist >> ch; // '{'
					  LoadVMatrixLine(ist, tempoldutofreezone.functorMethod, 3 * freezone.Size() - 1);
        
					  ist >> ch; // '{'
					  LoadVMatrixLine(ist, tempoldutofreezone.functorMethod, 3 * freezone.Size());
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
			  else if (string.Compare(buf, "freezone2") == 0)
			  {
			  int k;
			  int nfp;
        
			  nfp = 0;
			  ist >> ch;
        
			  DenseMatrix hm1 = new DenseMatrix(3, 50);
			  DenseMatrix hm2 = new DenseMatrix(50, 50);
			  DenseMatrix hm3 = new DenseMatrix(50, 50);
			  hm3 = 0;
        
			  while (ch == '{')
			  {
				  hm1 = 0;
				  nfp++;
				  LoadVMatrixLine(ist, hm1.functorMethod, 1);
        
				  for (i = 1; i <= points.Size(); i++)
				  {
				tfz.Elem(nfp, i) = hm1.Get(1, 3 * i - 2);
				  }
        
        
				  p.X() = p.Y() = p.Z() = 0;
				  for (i = 1; i <= points.Size(); i++)
				  {
				  p.X() += hm1.Get(1, 3 * i - 2) * points.Get(i).X();
				  p.Y() += hm1.Get(1, 3 * i - 2) * points.Get(i).Y();
				  p.Z() += hm1.Get(1, 3 * i - 2) * points.Get(i).Z();
				  }
				  freezone.Append(p);
				  freezonelimit.Append(p);
        
				  hm2 = 0;
				  for (i = 1; i <= 3 * noldp; i++)
				  {
				hm2.Elem(i, i) = 1;
				  }
				  for (i = 1; i <= 3 * noldp; i++)
				  {
				for (j = 1; j <= 3 * (points.Size() - noldp); j++)
				{
				  hm2.Elem(j + 3 * noldp, i) = tempoldutonewu.Get(j, i);
				}
				  }
        
				  for (i = 1; i <= 3; i++)
				  {
				for (j = 1; j <= 3 * noldp; j++)
				{
					double sum = 0;
					for (k = 1; k <= 3 * points.Size(); k++)
					{
					  sum += hm1.Get(i, k) * hm2.Get(k, j);
					}
        
					hm3.Elem(i + 3 * (nfp - 1), j) = sum;
				}
				  }
        
				  //	    (*testout) << "freepoint: " << p << endl;
        
				  while (ch != ';')
				  {
				ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  tfzl = tfz.functorMethod;
        
			  tempoldutofreezone.CopyFrom(hm3.functorMethod);
			  tempoldutofreezonelimit.CopyFrom(hm3.functorMethod);
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "freezonelimit") == 0)
			  {
			  int k;
			  int nfp;
			  nfp = 0;
			  ist >> ch;
        
			  DenseMatrix hm1 = new DenseMatrix(3, 50);
			  DenseMatrix hm2 = new DenseMatrix(50, 50);
			  DenseMatrix hm3 = new DenseMatrix(50, 50);
			  hm3 = 0;
        
			  while (ch == '{')
			  {
				  hm1 = 0;
				  nfp++;
				  LoadVMatrixLine(ist, hm1.functorMethod, 1);
        
				  for (i = 1; i <= points.Size(); i++)
				  {
				tfzl.Elem(nfp, i) = hm1.Get(1, 3 * i - 2);
				  }
        
        
				  p.X() = p.Y() = p.Z() = 0;
				  for (i = 1; i <= points.Size(); i++)
				  {
				  p.X() += hm1.Get(1, 3 * i - 2) * points.Get(i).X();
				  p.Y() += hm1.Get(1, 3 * i - 2) * points.Get(i).Y();
				  p.Z() += hm1.Get(1, 3 * i - 2) * points.Get(i).Z();
				  }
				  freezonelimit.Elem(nfp) = p;
        
				  hm2 = 0;
				  for (i = 1; i <= 3 * noldp; i++)
				  {
				hm2.Elem(i, i) = 1;
				  }
				  for (i = 1; i <= 3 * noldp; i++)
				  {
				for (j = 1; j <= 3 * (points.Size() - noldp); j++)
				{
				  hm2.Elem(j + 3 * noldp, i) = tempoldutonewu.Get(j, i);
				}
				  }
        
				  for (i = 1; i <= 3; i++)
				  {
				for (j = 1; j <= 3 * noldp; j++)
				{
					double sum = 0;
					for (k = 1; k <= 3 * points.Size(); k++)
					{
					  sum += hm1.Get(i, k) * hm2.Get(k, j);
					}
        
					hm3.Elem(i + 3 * (nfp - 1), j) = sum;
				}
				  }
        
				  //	    (*testout) << "freepoint: " << p << endl;
        
				  while (ch != ';')
				  {
				ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  tempoldutofreezonelimit.CopyFrom(hm3.functorMethod);
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "freeset") == 0)
			  {
			  freesets.Append(new Array<int>());
        
			  ist >> ch;
        
			  while (ch != ';')
			  {
				  ist.putback(ch);
				  ist >> i;
				  freesets.Last().Append(i);
				  ist >> ch;
			  }
			  }
        
			  else if (string.Compare(buf, "elements") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  elements.Append(new Element(ELEMENT_TYPE.TET));
        
				  //	      elements.Last().SetNP(1);
				  ist >> elements.Last().PNum(1);
				  ist >> ch; // ','
        
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(2);
				  ist >> elements.Last().PNum(2);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(3);
				  ist >> elements.Last().PNum(3);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(4);
				  elements.Last().SetType(ELEMENT_TYPE.TET);
				  ist >> elements.Last().PNum(4);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(5);
				  elements.Last().SetType(ELEMENT_TYPE.PYRAMID);
				  ist >> elements.Last().PNum(5);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(6);
				  elements.Last().SetType(ELEMENT_TYPE.PRISM);
				  ist >> elements.Last().PNum(6);
				  ist >> ch; // ','
				  }
        
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(6);
				  elements.Last().SetType(ELEMENT_TYPE.HEX);
				  ist >> elements.Last().PNum(7);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  //		  elements.Last().SetNP(6);
				  elements.Last().SetType(ELEMENT_TYPE.HEX);
				  ist >> elements.Last().PNum(8);
				  ist >> ch; // ','
				  }
        
				  /*
				  orientations.Append (fourint());
				  orientations.Last().i1 = elements.Last().PNum(1);
				  orientations.Last().i2 = elements.Last().PNum(2);
				  orientations.Last().i3 = elements.Last().PNum(3);
				  orientations.Last().i4 = elements.Last().PNum(4);
				  */
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "orientations") == 0)
        
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  //        fourint a = fourint();
				  orientations.Append(new fourint());
        
				  ist >> orientations.Last().i1;
				  ist >> ch; // ','
				  ist >> orientations.Last().i2;
				  ist >> ch; // ','
				  ist >> orientations.Last().i3;
				  ist >> ch; // ','
				  ist >> orientations.Last().i4;
				  ist >> ch; // ','
        
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
        
			  else if (string.Compare(buf, "endrule") != 0)
			  {
			  PrintSysError("Parser3d, unknown token ", buf);
			  }
		  } while (!ist.eof() && string.Compare(buf, "endrule") != 0);
        
        
		  //  (*testout) << endl;
		  //  (*testout) << Name() << endl;
		  //  (*testout) << "no1 = " << GetNO() << endl;
        
		  oldutonewu.SetSize(3 * (points.Size() - noldp), 3 * noldp);
		  oldutonewu = 0;
        
		  for (i = 1; i <= oldutonewu.Height(); i++)
		  {
			for (j = 1; j <= oldutonewu.Width(); j++)
			{
			  oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);
			}
		  }
        
        
		  /*
		    oldutofreezone = new SparseMatrixFlex (3 * freezone.Size(), 3 * noldp);
		    oldutofreezonelimit = new SparseMatrixFlex (3 * freezone.Size(), 3 * noldp);
		
		    oldutofreezone -> SetSymmetric(0);
		    oldutofreezonelimit -> SetSymmetric(0);
		    */
        
		  /*
		    oldutofreezone = new DenseMatrix (3 * freezone.Size(), 3 * noldp);
		    oldutofreezonelimit = new DenseMatrix (3 * freezone.Size(), 3 * noldp);
		  
		    for (i = 1; i <= oldutofreezone->Height(); i++)
		    for (j = 1; j <= oldutofreezone->Width(); j++)
		    //      if (j == 4 || j >= 7)
		    {
		    if (tempoldutofreezone.Elem(i, j))
		    (*oldutofreezone)(i, j) = tempoldutofreezone(i, j);
		    if (tempoldutofreezonelimit.Elem(i, j))
		    (*oldutofreezonelimit)(i, j) = tempoldutofreezonelimit(i, j);
		    }
		    */
        
        
        
        
		  oldutofreezone = new DenseMatrix(freezone.Size(), points.Size());
		  oldutofreezonelimit = new DenseMatrix(freezone.Size(), points.Size());
		  //  oldutofreezone = new SparseMatrixFlex (freezone.Size(), points.Size());
		  //  oldutofreezonelimit = new SparseMatrixFlex (freezone.Size(), points.Size());
        
		  for (i = 1; i <= freezone.Size(); i++)
		  {
			for (j = 1; j <= points.Size(); j++)
			{
			if (tfz.Elem(i, j) != 0)
			{
			  oldutofreezone.Elem(i, j) = tfz.Elem(i, j);
			}
			if (tfzl.Elem(i, j) != 0)
			{
			  oldutofreezonelimit.Elem(i, j) = tfzl.Elem(i, j);
			}
			}
		  }
        
		  /*
		  (*testout) << "Rule " << Name() << endl;
		  (*testout) << "oldutofreezone = " << (*oldutofreezone) << endl;
		  (*testout) << "oldutofreezonelimit = " << (*oldutofreezonelimit) << endl;
		  */
        
		  freezonepi.SetSize(freezone.Size());
		  for (i = 1; i <= freezonepi.Size(); i++)
		  {
			freezonepi.Elem(i) = 0;
		  }
		  for (i = 1; i <= freezone.Size(); i++)
		  {
			for (j = 1; j <= noldp; j++)
			{
			  if (Dist(freezone.Get(i), points.Get(j)) < 1e-8)
			  {
			freezonepi.Elem(i) = j;
			  }
			}
		  }
        
        
        
        
		  for (i = 1; i <= elements.Size(); i++)
		  {
			  if (elements.Elem(i).GetNP() == 4)
			  {
			  orientations.Append(new fourint());
			  orientations.Last().i1 = elements.Get(i).PNum(1);
			  orientations.Last().i2 = elements.Get(i).PNum(2);
			  orientations.Last().i3 = elements.Get(i).PNum(3);
			  orientations.Last().i4 = elements.Get(i).PNum(4);
			  }
			  if (elements.Elem(i).GetNP() == 5)
			  {
			  orientations.Append(new fourint());
			  orientations.Last().i1 = elements.Get(i).PNum(1);
			  orientations.Last().i2 = elements.Get(i).PNum(2);
			  orientations.Last().i3 = elements.Get(i).PNum(3);
			  orientations.Last().i4 = elements.Get(i).PNum(5);
        
			  orientations.Append(new fourint());
			  orientations.Last().i1 = elements.Get(i).PNum(1);
			  orientations.Last().i2 = elements.Get(i).PNum(3);
			  orientations.Last().i3 = elements.Get(i).PNum(4);
			  orientations.Last().i4 = elements.Get(i).PNum(5);
			  }
		  }
        
        
        
		  if (freesets.Size() == 0)
		  {
			  freesets.Append(new Array<int>());
			  for (i = 1; i <= freezone.Size(); i++)
			  {
			freesets.Elem(1).Append(i);
			  }
		  }
        
        
		  //  testout << "Freezone: " << endl;
        
		  //  for (i = 1; i <= freezone.Size(); i++)
		  //    (*testout) << "freepoint: " << freezone.Get(i) << endl;
		  Vector vp = new Vector(points.Size());
		  Vector vfp = new Vector(freezone.Size());
        
        
		  if (quality < 100)
		  {
			  for (int i = 1; i <= 3; i++)
			  {
			  for (int j = 1; j <= points.Size(); j++)
			  {
				vp(j - 1) = points.Get(j).X(i);
			  }
			  oldutofreezone.Mult(vp, vfp);
			  for (int j = 1; j <= freezone.Size(); j++)
			  {
				freezone.Elem(j).X(i) = vfp(j - 1);
			  }
			  }
			  //      for (i = 1; i <= freezone.Size(); i++)
			  //	(*testout) << "freepoint: " << freezone.Get(i) << endl;
		  }
        
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  freefaces.Append(new Array<threeint>());
        
			  Array<int> freeset = *freesets.Elem(fs);
			  Array<threeint> freesetfaces = *freefaces.Last();
        
			  for (ii1 = 1; ii1 <= freeset.Size(); ii1++)
			  {
			for (ii2 = 1; ii2 <= freeset.Size(); ii2++)
			{
			  for (ii3 = 1; ii3 <= freeset.Size(); ii3++)
			  {
				if (ii1 < ii2 && ii1 < ii3 && ii2 != ii3)
				{
				i1 = freeset.Get(ii1);
				i2 = freeset.Get(ii2);
				i3 = freeset.Get(ii3);
        
				Vec3d v1 = new Vec3d();
				Vec3d v2 = new Vec3d();
				Vec3d n = new Vec3d();
        
				v1 = freezone.Get(i3) - freezone.Get(i1);
				v2 = freezone.Get(i2) - freezone.Get(i1);
				n = Cross(new netgen.Vec3d(v1), new netgen.Vec3d(v2));
				n /= n.Length();
				//		(*testout) << "i1,2,3 = " << i1 << ", " << i2 << ", " << i3 << endl;
				//		(*testout) << "v1 = " << v1 << " v2 = " << v2 << " n = " << n << endl;
				ok = 1;
				for (ii = 1; ii <= freeset.Size(); ii++)
				{
					i = freeset.Get(ii);
					//		    (*testout) << "i = " << i << endl;
					if (i != i1 && i != i2 && i != i3)
					{
					  if ((freezone.Get(i) - freezone.Get(i1)) * n < 0)
					  {
						  ok = 0;
					  }
					}
				}
        
				if (ok)
				{
					freesetfaces.Append(new threeint());
					freesetfaces.Last().i1 = i1;
					freesetfaces.Last().i2 = i2;
					freesetfaces.Last().i3 = i3;
				}
				}
			  }
			}
			  }
		  }
        
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  freefaceinequ.Append(new DenseMatrix(freefaces.Get(fs).Size(), 4));
		  }
        
        
		  {
			int minn;
			//    Array<int> pnearness (noldp);
			pnearness.SetSize(noldp);
        
			for (i = 1; i <= pnearness.Size(); i++)
			{
			  pnearness.Elem(i) = INT_MAX / 10;
			}
        
			for (j = 1; j <= GetNP(1); j++)
			{
			  pnearness.Elem(GetPointNr(1, j)) = 0;
			}
        
			do
			{
			ok = 1;
        
			for (i = 1; i <= noldf; i++)
			{
				minn = INT_MAX / 10;
				for (j = 1; j <= GetNP(i); j++)
				{
				  minn = min2(minn, pnearness.Get(GetPointNr(i, j)));
				}
        
				for (j = 1; j <= GetNP(i); j++)
				{
				  if (pnearness.Get(GetPointNr(i, j)) > minn + 1)
				  {
				  ok = 0;
				  pnearness.Elem(GetPointNr(i, j)) = minn + 1;
				  }
				}
			}
        
			for (i = 1; i <= edges.Size(); i++)
			{
				int pi1 = edges.Get(i).i1;
				int pi2 = edges.Get(i).i2;
        
				if (pnearness.Get(pi1) > pnearness.Get(pi2) + 1)
				{
				ok = 0;
				pnearness.Elem(pi1) = pnearness.Get(pi2) + 1;
				}
				if (pnearness.Get(pi2) > pnearness.Get(pi1) + 1)
				{
				ok = 0;
				pnearness.Elem(pi2) = pnearness.Get(pi1) + 1;
				}
			}
        
        
			for (i = 1; i <= elements.Size(); i++)
			{
			  if (elements.Get(i).GetNP() == 6) // prism rule
			  {
				  for (j = 1; j <= 3; j++)
				  {
				  int pi1 = elements.Get(i).PNum(j);
				  int pi2 = elements.Get(i).PNum(j + 3);
        
				  if (pnearness.Get(pi1) > pnearness.Get(pi2) + 1)
				  {
					  ok = 0;
					  pnearness.Elem(pi1) = pnearness.Get(pi2) + 1;
				  }
				  if (pnearness.Get(pi2) > pnearness.Get(pi1) + 1)
				  {
					  ok = 0;
					  pnearness.Elem(pi2) = pnearness.Get(pi1) + 1;
				  }
				  }
			  }
			}
			} while (!ok);
        
			maxpnearness = 0;
			for (i = 1; i <= pnearness.Size(); i++)
			{
			  maxpnearness = max2(maxpnearness, pnearness.Get(i));
			}
        
        
			fnearness.SetSize(noldf);
        
			for (i = 1; i <= noldf; i++)
			{
			fnearness.Elem(i) = 0;
			for (j = 1; j <= GetNP(i); j++)
			{
			  fnearness.Elem(i) += pnearness.Get(GetPointNr(i, j));
			}
			}
        
			// (*testout) << "rule " << name << ", pnear = " << pnearness << endl;
		  }
        
        
		  //Table of edges:
		  for (fs = 1; fs <= freesets.Size(); fs++)
		  {
			  freeedges.Append(new Array<twoint>());
        
			  //      Array<int> & freeset = *freesets.Get(fs);
			  Array<twoint> freesetedges = *freeedges.Last();
			  Array<threeint> freesetfaces = *freefaces.Get(fs);
			  int k;
			  int l;
			  int ind;
        
			  for (k = 1; k <= freesetfaces.Size(); k++)
			  {
				  // threeint tr = freesetfaces.Get(k);
        
			  for (l = k + 1; l <= freesetfaces.Size(); l++)
			  {
				  ind = NeighbourTrianglePoint(freesetfaces.Get(k), freesetfaces.Get(l));
				  if (ind == 0)
				  {
					  continue;
				  }
        
				  INDEX_3 f1 = new INDEX_3(freesetfaces.Get(k).i1, freesetfaces.Get(k).i2, freesetfaces.Get(k).i3);
				  INDEX_3 f2 = new INDEX_3(freesetfaces.Get(l).i1, freesetfaces.Get(l).i2, freesetfaces.Get(l).i3);
				  INDEX_2 ed = new INDEX_2(0, 0);
				  for (int f11 = 1; f11 <= 3; f11++)
				  {
				for (int f12 = 1; f12 <= 3; f12++)
				{
				  if (f11 != f12)
				  {
					for (int f21 = 1; f21 <= 3; f21++)
					{
					  for (int f22 = 1; f22 <= 3; f22++)
					  {
					if (f1.I(f11) == f2.I(f21) && f1.I(f12) == f2.I(f22))
					{
					  ed.I(1) = f1.I(f11);
					  ed.I(2) = f1.I(f12);
					}
					  }
					}
				  }
				}
				  }
				  //	      (*testout) << "ed = " << ed.I(1) << "-" << ed.I(2) << endl;
				  //	      (*testout) << "ind = " << ind << " ed = " << ed << endl;
				  for (int eli = 1; eli <= GetNOldF(); eli++)
				  {
				  if (GetNP(eli) == 4)
				  {
					  for (int elr = 1; elr <= 4; elr++)
					  {
					  if (GetPointNrMod(eli, elr) == ed.I(1) && GetPointNrMod(eli, elr + 2) == ed.I(2))
					  {
						  /*
						  (*testout) << "ed is diagonal of rectangle" << endl;
						  (*testout) << "ed = " << ed.I(1) << "-" << ed.I(2) << endl;
						  (*testout) << "ind = " << ind << endl;
						  */
						  ind = 0;
					  }
        
					  }
				  }
				  }
        
				  if (ind != 0)
				  {
				  /*
				  (*testout) << "new edge from face " << k 
						 << " = (" << freesetfaces.Get(k).i1 
						 << ", " << freesetfaces.Get(k).i2 
						 << ", " << freesetfaces.Get(k).i3
						 << "), point " << ind << endl;
						 */
				  freesetedges.Append(new twoint(k, ind));
				  }
			  }
			  }
		  }
        
		}
	}
}