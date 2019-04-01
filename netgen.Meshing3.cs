namespace netgen
{

	public class Meshing3
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Delaunay(Mesh mesh, int domainnr, MeshingParameters mp)
		  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static Timer t("Meshing3::Delaunay");
			RegionTimer reg = new RegionTimer(Delaunay_t);
        
			int np;
			int ne;
        
			PrintMessage(1, "Delaunay meshing");
			PrintMessage(3, "number of points: ", mesh.GetNP());
			PushStatus("Delaunay meshing");
        
        
			Array<DelaunayTet> tempels = new Array<DelaunayTet>();
			Point3d pmin = new Point3d();
			Point3d pmax = new Point3d();
        
			DelaunayTet startel = new DelaunayTet();
        
			int oldnp = mesh.GetNP();
			if (mp.blockfill != 0)
			{
			BlockFillLocalH(mesh, mp);
			PrintMessage(3, "number of points: ", mesh.GetNP());
			}
        
			np = mesh.GetNP();
        
			Delaunay1(mesh, mp, adfront, tempels, oldnp, startel, ref pmin, ref pmax);
        
			{
			  // improve delaunay - mesh by swapping !!!!
        
			  Mesh tempmesh = new Mesh();
        
			  foreach (var meshpoint in mesh.Points())
			  {
				tempmesh.AddPoint(meshpoint);
			  }
        
			  foreach (var tempel in tempels)
			  {
			  Element el = new Element(4);
			  for (int j = 0; j < 4; j++)
			  {
					el[j] = tempel[j];
			  }
        
			  el.SetIndex(1);
        
			  const Point < 3> & lp1 = new mesh.Point(el[0]);
			  const Point < 3> & lp2 = new mesh.Point(el[1]);
			  const Point < 3> & lp3 = new mesh.Point(el[2]);
			  const Point < 3> & lp4 = new mesh.Point(el[3]);
			  Vec < 3> v1 = lp2 - lp1;
			  Vec < 3> v2 = lp3 - lp1;
			  Vec < 3> v3 = lp4 - lp1;
        
			  Vec < 3> n = Cross(v1, v2);
			  double vol = n * v3;
			  if (vol > 0)
			  {
				  swap(el[2], el[3]);
			  }
        
			  tempmesh.AddVolumeElement(el);
			  }
        
        
			  MeshQuality3d(tempmesh);
        
			  tempmesh.AddFaceDescriptor(new FaceDescriptor(1, 1, 0, 0));
			  tempmesh.AddFaceDescriptor(new FaceDescriptor(2, 1, 0, 0));
        
        
        
			  for (int i = 1; i <= mesh.GetNOpenElements(); i++)
			  {
			  Element2d sel = mesh.OpenElement(i);
			  sel.SetIndex(1);
			  tempmesh.AddSurfaceElement(sel);
			  swap(sel[1], sel[2]);
			  tempmesh.AddSurfaceElement(sel);
			  }
        
        
			  for (int i = 1; i <= 4; i++)
			  {
			  Element2d self = new Element2d(ELEMENT_TYPE.TRIG);
			  self.SetIndex(1);
			  startel.GetFace(i - 1, self);
			  tempmesh.AddSurfaceElement(self);
			  }
        
        
			  //  for (i = mesh.GetNP() - 3; i <= mesh.GetNP(); i++)
			  //    tempmesh.AddLockedPoint (i);
			  foreach (var pi in tempmesh.Points().Range())
			  {
				tempmesh.AddLockedPoint(pi);
			  }
        
			  //  tempmesh.PrintMemInfo(cout);
			  // tempmesh.Save ("tempmesh.vol");
        
			  for (int i = 1; i <= 4; i++)
			  {
			  tempmesh.FindOpenElements();
        
			  PrintMessage(5, "Num open: ", tempmesh.GetNOpenElements());
			  tempmesh.CalcSurfacesOfNode();
        
			  tempmesh.FreeOpenElementsEnvironment(1);
        
			  MeshOptimize3d meshopt = new MeshOptimize3d(mp);
			  // tempmesh.CalcSurfacesOfNode();
				  meshopt.SwapImprove(tempmesh, OPTIMIZEGOAL.OPT_CONFORM);
			  }
        
			  MeshQuality3d(tempmesh);
        
			  tempels.SetSize(0);
			  foreach (var el in tempmesh.VolumeElements())
			  {
				tempels.Append(el);
			  }
			}
        
        
        
			// remove degenerated
        
			BitArray badnode = new BitArray(mesh.GetNP());
			badnode.Clear();
			int ndeg = 0;
			for (int i = 1; i <= tempels.Size(); i++)
			{
			Element el = new Element(4);
			for (int j = 0; j < 4; j++)
			{
			  el[j] = tempels.Elem(i)[j];
			}
			//      Element & el = tempels.Elem(i);
			Point3d lp1 = new mesh.Point(el[0]);
			Point3d lp2 = new mesh.Point(el[1]);
			Point3d lp3 = new mesh.Point(el[2]);
			Point3d lp4 = new mesh.Point(el[3]);
			Vec3d v1 = new Vec3d(lp1, lp2);
			Vec3d v2 = new Vec3d(lp1, lp3);
			Vec3d v3 = new Vec3d(lp1, lp4);
			Vec3d n = Cross(new netgen.Vec3d(v1), new netgen.Vec3d(v2));
			double vol = n * v3;
        
			double h = v1.Length() + v2.Length() + v3.Length();
			if (ngsimd.GlobalMembers.fabs(vol) < 1e-8 * (h * h * h) && (el[0] <= np && el[1] <= np && el[2] <= np && el[3] <= np)) // old: 1e-12
			{
				badnode.Set(new netgen.Element(el[0]));
				badnode.Set(new netgen.Element(el[1]));
				badnode.Set(new netgen.Element(el[2]));
				badnode.Set(new netgen.Element(el[3]));
				ndeg++;
				(*testout) << "vol = " << vol << " h = " << h << "\n";
			}
        
			if (vol > 0)
			{
			  Swap(ref el[2], ref el[3]);
			}
			}
        
			ne = tempels.Size();
			for (int i = ne; i >= 1; i--)
			{
			DelaunayTet el = tempels.Get(i);
			if (badnode.Test(new netgen.DelaunayTet(el[0])) || badnode.Test(new netgen.DelaunayTet(el[1])) || badnode.Test(new netgen.DelaunayTet(el[2])) || badnode.Test(new netgen.DelaunayTet(el[3])))
			{
			  tempels.DeleteElement(i);
			}
			}
        
        
			PrintMessage(3, ndeg, " degenerated elements removed");
        
			// find surface triangles which are no face of any tet
        
			INDEX_3_HASHTABLE<int> openeltab = new INDEX_3_HASHTABLE<int>(mesh.GetNOpenElements() + 3);
			Array<int> openels = new Array<int>();
			for (int i = 1; i <= mesh.GetNOpenElements(); i++)
			{
			Element2d tri = mesh.OpenElement(i);
			INDEX_3 i3 = new INDEX_3(tri[0], tri[1], tri[2]);
			i3.Sort();
			openeltab.Set(i3, i);
			}
        
			for (int i = 1; i <= tempels.Size(); i++)
			{
			for (int j = 0; j < 4; j++)
			{
				INDEX_3 i3 = tempels.Get(i).GetFace(j);
				i3.Sort();
				if (openeltab.Used(i3))
				{
				  openeltab.Set(i3, 0);
				}
			}
			}
        
			// and store them in openels
			for (int i = 1; i <= openeltab.GetNBags(); i++)
			{
			  for (int j = 1; j <= openeltab.GetBagSize(i); j++)
			  {
			  INDEX_3 i3 = new INDEX_3();
			  int fnr;
			  openeltab.GetData(i, j, i3, fnr);
			  if (fnr != 0)
			  {
				openels.Append(fnr);
			  }
			  }
			}
        
        
        
        
        
			// find open triangle with close edge (from halfening of surface squares)
        
			INDEX_2_HASHTABLE<INDEX_2> twotrias = new INDEX_2_HASHTABLE<INDEX_2>(mesh.GetNOpenElements() + 5);
			//  for (i = 1; i <= mesh.GetNOpenElements(); i++)
			for (int ii = 1; ii <= openels.Size(); ii++)
			{
			int i = openels.Get(ii);
			Element2d el = mesh.OpenElement(i);
			for (int j = 1; j <= 3; j++)
			{
				INDEX_2 hi2 = new INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
				hi2.Sort();
				if (twotrias.Used(hi2))
				{
				INDEX_2 hi3 = new INDEX_2();
				hi3 = twotrias.Get(hi2);
				hi3.I2() = el.PNumMod(j + 2);
				twotrias.Set(hi2, hi3);
				}
				else
				{
				INDEX_2 hi3 = new INDEX_2(el.PNumMod(j + 2), 0);
				twotrias.Set(hi2, hi3);
				}
			}
			}
        
			INDEX_2_HASHTABLE<int> tetedges = new INDEX_2_HASHTABLE<int>(tempels.Size() + 5);
			for (int i = 1; i <= tempels.Size(); i++)
			{
			DelaunayTet el = tempels.Get(i);
			INDEX_2 i2 = new INDEX_2();
			for (int j = 1; j <= 6; j++)
			{
				switch (j)
				{
				  case 1:
					  i2.I1() = el[0];
					  i2.I2() = el[1];
					  break;
				  case 2:
					  i2.I1() = el[0];
					  i2.I2() = el[2];
					  break;
				  case 3:
					  i2.I1() = el[0];
					  i2.I2() = el[3];
					  break;
				  case 4:
					  i2.I1() = el[1];
					  i2.I2() = el[2];
					  break;
				  case 5:
					  i2.I1() = el[1];
					  i2.I2() = el[3];
					  break;
				  case 6:
					  i2.I1() = el[2];
					  i2.I2() = el[3];
					  break;
				  default:
					  i2.I1() = i2.I2() = 0;
					  break;
				}
				i2.Sort();
				tetedges.Set(i2, 1);
			}
			}
			//  cout << "tetedges:";
			//  tetedges.PrintMemInfo (cout);
        
        
			for (INDEX_2_HASHTABLE<INDEX_2>.Iterator it = twotrias.Begin(); it != twotrias.End(); it++)
			{
			INDEX_2 hi2 = new INDEX_2();
			INDEX_2 hi3 = new INDEX_2();
			twotrias.GetData(it, ref hi2, ref hi3);
			hi3.Sort();
			if (tetedges.Used(hi3))
			{
				Point3d p1 = new mesh.Point(new PointIndex(hi2.I1()));
				Point3d p2 = new mesh.Point(new PointIndex(hi2.I2()));
				Point3d p3 = new mesh.Point(new PointIndex(hi3.I1()));
				Point3d p4 = new mesh.Point(new PointIndex(hi3.I2()));
				Vec3d v1 = new Vec3d(p1, p2);
				Vec3d v2 = new Vec3d(p1, p3);
				Vec3d v3 = new Vec3d(p1, p4);
				Vec3d n = Cross(new netgen.Vec3d(v1), new netgen.Vec3d(v2));
				double vol = n * v3;
        
				double h = v1.Length() + v2.Length() + v3.Length();
				if (ngsimd.GlobalMembers.fabs(vol) < 1e-4 * (h * h * h)) // old: 1e-12
				{
				badnode.Set(hi3.I1());
				badnode.Set(hi3.I2());
				}
			}
			}
        
			/*
			for (i = 1; i <= twotrias.GetNBags(); i++)
			  for (j = 1; j <= twotrias.GetBagSize (i); j++)
			{
			  INDEX_2 hi2, hi3;
			  twotrias.GetData (i, j, hi2, hi3);
			  hi3.Sort();
			  if (tetedges.Used (hi3))
				{
				  const Point3d & p1 = mesh.Point (hi2.I1());
				  const Point3d & p2 = mesh.Point (hi2.I2());
				  const Point3d & p3 = mesh.Point (hi3.I1());
				  const Point3d & p4 = mesh.Point (hi3.I2());
				  Vec3d v1(p1, p2);
				  Vec3d v2(p1, p3);
				  Vec3d v3(p1, p4);
				  Vec3d n = Cross (v1, v2);
				  double vol = n * v3;
		
				  double h = v1.Length() + v2.Length() + v3.Length();
				  if (fabs (vol) < 1e-4 * (h * h * h))   // old: 1e-12
				{
				  badnode.Set(hi3.I1());	
				  badnode.Set(hi3.I2());	
				}
				}
			}
			*/
        
			ne = tempels.Size();
			for (int i = ne; i >= 1; i--)
			{
			DelaunayTet el = tempels.Get(i);
			if (badnode.Test(new netgen.DelaunayTet(el[0])) || badnode.Test(new netgen.DelaunayTet(el[1])) || badnode.Test(new netgen.DelaunayTet(el[2])) || badnode.Test(new netgen.DelaunayTet(el[3])))
			{
			  tempels.DeleteElement(i);
			}
			}
        
        
        
        
			// find intersecting:
			PrintMessage(3, "Remove intersecting");
			if (openels.Size())
			{
			BoxTree < 3> setree(pmin, pmax);
        
			/*      
				cout << "open elements in search tree: " << openels.Size() << endl;
				cout << "pmin, pmax = " << pmin << " - " << pmax << endl;
			*/
        
			for (int i = 1; i <= openels.Size(); i++)
			{
				int fnr;
				fnr = openels.Get(i);
				if (fnr != 0)
				{
				Element2d tri = mesh.OpenElement(fnr);
        
				Point3d ltpmin = new Point3d(new mesh.Point(tri[0]));
				Point3d ltpmax = new Point3d(ltpmin);
        
				for (int k = 2; k <= 3; k++)
				{
					ltpmin.SetToMin(new mesh.Point(tri.PNum(k)));
					ltpmax.SetToMax(new mesh.Point(tri.PNum(k)));
				}
				setree.Insert(ltpmin, ltpmax, fnr);
				}
			}
        
			Array<int> neartrias = new Array<int>();
			for (int i = 1; i <= tempels.Size(); i++)
			{
				const Point < 3> *pp[4];
				int[] tetpi = new int[4];
				DelaunayTet el = tempels.Elem(i);
        
				int intersect = 0;
        
				for (int j = 0; j < 4; j++)
				{
				pp[j] = &new mesh.Point(el[j]);
				tetpi[j] = el[j];
				}
        
				Point3d tetpmin = new Point3d(*pp[0]);
				Point3d tetpmax = new Point3d(tetpmin);
				for (int j = 1; j < 4; j++)
				{
				tetpmin.SetToMin(pp[j]);
				tetpmax.SetToMax(pp[j]);
				}
				tetpmin.CopyFrom(tetpmin + 0.01 * (tetpmin - tetpmax));
				tetpmax.CopyFrom(tetpmax + 0.01 * (tetpmax - tetpmin));
        
				setree.GetIntersecting(tetpmin, tetpmax, neartrias);
        
        
				//      for (j = 1; j <= mesh.GetNSE(); j++)
				//	{
				for (int jj = 1; jj <= neartrias.Size(); jj++)
				{
				int j = neartrias.Get(jj);
        
				Element2d tri = mesh.OpenElement(j);
				const Point < 3> *tripp[3];
				int[] tripi = new int[3];
        
				for (int k = 1; k <= 3; k++)
				{
					tripp[k - 1] = &new mesh.Point(tri.PNum(k));
					tripi[k - 1] = tri.PNum(k);
				}
        
				if (IntersectTetTriangle(pp[0], tripp[0], tetpi, tripi))
				{
					/*
					int il1, il2;
					(*testout) << "intersect !" << endl;
					(*testout) << "triind: ";
					for (il1 = 0; il1 < 3; il1++)
					  (*testout) << " " << tripi[il1];
					(*testout) << endl;
					(*testout) << "tetind: ";
					for (il2 = 0; il2 < 4; il2++)
					  (*testout) << " " << tetpi[il2];
					(*testout) << endl;
		
					(*testout) << "trip: ";
					for (il1 = 0; il1 < 3; il1++)
					  (*testout) << " " << *tripp[il1];
					(*testout) << endl;
					(*testout) << "tetp: ";
					for (il2 = 0; il2 < 4; il2++)
					  (*testout) << " " << *pp[il2];
					(*testout) << endl;
					*/
        
        
					intersect = 1;
					break;
				}
				}
        
        
				if (intersect != 0)
				{
				tempels.DeleteElement(i);
				i--;
				}
			}
			}
        
        
        
        
			PrintMessage(3, "Remove outer");
        
			// find connected tets (with no face between, and no hole due
			// to removed intersecting tets.
			//  INDEX_3_HASHTABLE<INDEX_2> innerfaces(np);
        
        
			INDEX_3_HASHTABLE<int> boundaryfaces = new INDEX_3_HASHTABLE<int>(mesh.GetNOpenElements() / 3 + 1);
			/*
			for (int i = 1; i <= mesh.GetNOpenElements(); i++)
			  {
			const Element2d & tri = mesh.OpenElement(i);
			INDEX_3 i3 (tri[0], tri[1], tri[2]);
			i3.Sort();
			boundaryfaces.PrepareSet (i3);
			  }
			*/
			foreach (Element2d tri in mesh.OpenElements())
			{
			INDEX_3 i3 = new INDEX_3(tri[0], tri[1], tri[2]);
			i3.Sort();
			boundaryfaces.PrepareSet(i3);
			}
			boundaryfaces.AllocateElements();
			for (int i = 1; i <= mesh.GetNOpenElements(); i++)
			{
			Element2d tri = mesh.OpenElement(i);
			INDEX_3 i3 = new INDEX_3(tri[0], tri[1], tri[2]);
			i3.Sort();
			boundaryfaces.Set(i3, 1);
			}
        
			/*
			for (int i = 0; i < tempels.Size(); i++)
			  for (int j = 0; j < 4; j++)
			tempels[i].NB(j) = 0;
			*/
			foreach (var el in tempels)
			{
			  for (int j = 0; j < 4; j++)
			  {
			el.NB(j) = 0;
			  }
			}
        
			TABLE<int,PointIndex.BASE> elsonpoint = new TABLE<int,PointIndex.BASE>(mesh.GetNP());
			/*
			for (int i = 0; i < tempels.Size(); i++)
			  {
			const DelaunayTet & el = tempels[i];
			*/
			foreach (DelaunayTet el in tempels)
			{
			INDEX_4 i4 = new INDEX_4(el[0], el[1], el[2], el[3]);
			i4.Sort();
			elsonpoint.IncSizePrepare(i4.I1());
			elsonpoint.IncSizePrepare(i4.I2());
			}
        
			elsonpoint.AllocateElementsOneBlock();
        
			for (int i = 0; i < tempels.Size(); i++)
			{
			DelaunayTet el = tempels[i];
			INDEX_4 i4 = new INDEX_4(el[0], el[1], el[2], el[3]);
			i4.Sort();
			elsonpoint.Add(i4.I1(), i + 1);
			elsonpoint.Add(i4.I2(), i + 1);
			}
        
			//  cout << "elsonpoint mem: ";
			//  elsonpoint.PrintMemInfo(cout);
        
			INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht = new INDEX_3_CLOSED_HASHTABLE<INDEX_2>(100);
        
			Element2d hel = new Element2d(ELEMENT_TYPE.TRIG);
			// for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
			foreach (PointIndex pi in mesh.Points().Range())
			{
			faceht.SetSize(4 * elsonpoint[pi].Size());
			for (int ii = 0; ii < elsonpoint[pi].Size(); ii++)
			{
				int i = elsonpoint[pi][ii];
				DelaunayTet el = tempels.Get(i);
        
				for (int j = 1; j <= 4; j++)
				{
				el.GetFace(j - 1, hel);
				hel.Invert();
				hel.NormalizeNumbering();
        
				if (hel[0] == pi)
				{
					INDEX_3 i3 = new INDEX_3(hel[0], hel[1], hel[2]);
        
					if (!boundaryfaces.Used(i3))
					{
					if (faceht.Used(i3))
					{
						INDEX_2 i2 = faceht.Get(i3);
        
						tempels.Elem(i).NB(j - 1) = i2.I1();
						tempels.Elem(i2.I1()).NB(i2.I2() - 1) = i;
					}
					else
					{
						hel.Invert();
						hel.NormalizeNumbering();
						INDEX_3 i3i = new INDEX_3(hel[0], hel[1], hel[2]);
						INDEX_2 i2 = new INDEX_2(i, j);
						faceht.Set(i3i, i2);
					}
					}
				}
				}
			}
			}
        
			/*
			  for (i = 1; i <= tempels.Size(); i++)
			  {
			  const DelaunayTet & el = tempels.Get(i);
			  for (j = 1; j <= 4; j++)
			  {
			  INDEX_3 i3;
			  Element2d face;
			  el.GetFace1 (j, face);
			  for (int kk = 1; kk <= 3; kk++)
			  i3.I(kk) = face.PNum(kk);
		
			  i3.Sort();
			  if (!boundaryfaces.Used (i3))
			  {
			  if (innerfaces.Used(i3))
			  {
			  INDEX_2 i2;
			  i2 = innerfaces.Get(i3);
			  i2.I2() = i;
			  innerfaces.Set (i3, i2);
			  }
			  else
			  {
			  INDEX_2 i2;
			  i2.I1() = i;
			  i2.I2() = 0;
			  innerfaces.Set (i3, i2);
			  }
			  }
			  }
			  }
			*/
        
			/*
			  (*testout) << "nb elements:" << endl;
			  for (i = 1; i <= tempels.Size(); i++)
			  {
			  (*testout) << i << " ";
			  for (j = 1; j <= 4; j++)
			  (*testout) << tempels.Get(i).NB1(j) << " ";
			  (*testout) << endl;
			  }
		
			  (*testout) << "pairs:" << endl;
			  for (i = 1; i <= innerfaces.GetNBags(); i++)
			  for (j = 1; j <= innerfaces.GetBagSize(i); j++)
			  {
			  INDEX_3 i3;
			  INDEX_2 i2;
			  innerfaces.GetData (i, j, i3, i2);
			  (*testout) << i2 << endl;
			  }
			*/
        
        
        
        
        
        
        
			/*
			  cout << "innerfaces: ";
			  innerfaces.PrintMemInfo (cout);
			*/
        
			//  cout << "boundaryfaces: ";
			//  boundaryfaces.PrintMemInfo (cout);
        
        
			PrintMessage(5, "tables filled");
        
        
			ne = tempels.Size();
			BitArray inner = new BitArray(ne);
			BitArray outer = new BitArray(ne);
			inner.Clear();
			outer.Clear();
			Array<int> elstack = new Array<int>();
        
			/*
			  int starti = 0;
			  for (i = 1; i <= ne; i++)
			  {
			  const Element & el = tempels.Get(i);
			  for (j = 1; j <= 4; j++)
			  for (k = 1; k <= 4; k++)
			  if (el.PNum(j) == startel.PNum(k))
			  {
			  outer.Set(i);
			  starti = i;
			  }
			  }
			*/
        
			while (true)
			{
			int inside;
			bool done = true;
        
			int i;
			for (i = 1; i <= ne; i++)
			{
			  if (!inner.Test(i) && !outer.Test(i))
			  {
				  done = false;
				  break;
			  }
			}
        
			if (done)
			{
				break;
			}
        
			DelaunayTet el = tempels.Get(i);
			Point3d p1 = new mesh.Point(el[0]);
			Point3d p2 = new mesh.Point(el[1]);
			Point3d p3 = new mesh.Point(el[2]);
			Point3d p4 = new mesh.Point(el[3]);
        
			Point3d ci = Center(p1, p2, p3, p4);
        
			inside = adfront.Inside(ci);
        
			/*
			  cout << "startel: " << i << endl;
			  cout << "inside = " << inside << endl;
			  cout << "ins2 = " << adfront->Inside (Center (ci, p1)) << endl;
			  cout << "ins3 = " << adfront->Inside (Center (ci, p2)) << endl;
			*/
        
			elstack.SetSize(0);
			elstack.Append(i);
        
			while (elstack.Size())
			{
				int ei = elstack.Last();
				elstack.DeleteLast();
        
				if (!inner.Test(ei) && !outer.Test(ei))
				{
				if (inside != 0)
				{
				  inner.Set(ei);
				}
				else
				{
				  outer.Set(ei);
				}
        
        
				for (int j = 1; j <= 4; j++)
				{
					INDEX_3 i3 = tempels.Get(ei).GetFace(j - 1);
					/*
					Element2d face;
					tempels.Get(ei).GetFace(j, face);
					for (int kk = 1; kk <= 3; kk++)
					  i3.I(kk) = face.PNum(kk);
					*/
					i3.Sort();
        
        
					if (tempels.Get(ei).NB(j - 1))
					{
					  elstack.Append(tempels.Get(ei).NB(j - 1));
					}
        
					/*
					  if (innerfaces.Used(i3))
					  {
					  INDEX_2 i2 = innerfaces.Get(i3);
					  int other = i2.I1() + i2.I2() - ei;
		
					  if (other != tempels.Get(ei).NB1(j))
					  cerr << "different1 !!" << endl;
		
					  if (other)
					  {
					  elstack.Append (other);
					  }
					  }
					  else
					  if (tempels.Get(ei).NB1(j))
					  cerr << "different2 !!" << endl;
					*/
        
				}
				}
			}
			}
        
        
        
			// check outer elements
			if (debugparam.slowchecks)
			{
			for (int i = 1; i <= ne; i++)
			{
				DelaunayTet el = tempels.Get(i);
				Point3d p1 = new mesh.Point(el[0]);
				Point3d p2 = new mesh.Point(el[1]);
				Point3d p3 = new mesh.Point(el[2]);
				Point3d p4 = new mesh.Point(el[3]);
        
				Point3d ci = Center(p1, p2, p3, p4);
        
				//       if (adfront->Inside (ci) != adfront->Inside (Center (ci, p1)))
				// 	cout << "ERROR: outer test unclear !!!" << endl;	
        
				if (inner.Test(i) != adfront.Inside(ci))
				{
				/*
				  cout << "ERROR: outer test wrong !!!" 
				  << "inner = " << int(inner.Test(i))
				  << "outer = " << int(outer.Test(i))
				  << endl;
		
				  cout << "Vol = " << Determinant(Vec3d(p1, p2),
				  Vec3d(p1, p3),
				  Vec3d(p1, p4)) << endl;
		
				*/	      
				for (int j = 1; j <= 4; j++)
				{
					Point3d hp = new Point3d();
					switch (j)
					{
					  case 1:
						  hp = Center(ci, p1);
						  break;
					  case 2:
						  hp = Center(ci, p2);
						  break;
					  case 3:
						  hp = Center(ci, p3);
						  break;
					  case 4:
						  hp = Center(ci, p4);
						  break;
					}
					//		  cout << "inside(" << hp << ") = " << adfront->Inside(hp) << endl;
				}
        
				}
        
				if (adfront.Inside(ci))
				{
				  outer.Clear(i);
				}
				else
				{
				  outer.Set(i);
				}
			}
			}
        
        
			/*
		
			// find bug in innerfaces
		
			tempmesh.DeleteVolumeElements();
		
			for (i = 1; i <= innerfaces.GetNBags(); i++)
			for (j = 1; j <= innerfaces.GetBagSize(i); j++)
			{
			INDEX_3 i3;
			INDEX_2 i2;
			innerfaces.GetData (i, j, i3, i2);
			if (i2.I2())
			{
			if (outer.Test(i2.I1()) != outer.Test(i2.I2()))
			{
			tempmesh.AddVolumeElement (tempels.Get(i2.I1()));
			tempmesh.AddVolumeElement (tempels.Get(i2.I2()));
			cerr << "outer flag different for connected els" << endl;
			}
			}
			}
		
		
			cout << "Check intersectiong once more" << endl;
		
			for (i = 1; i <= openels.Size(); i++)
			{
			tempmesh.SurfaceElement(2*openels.Get(i)).SetIndex(2);
			tempmesh.SurfaceElement(2*openels.Get(i)-1).SetIndex(2);
			}
		
			//  for (i = 1; i <= tempmesh.GetNE(); i++)
			//    for (j = 1; j <= tempmesh.GetNSE(); j++)
			i = 6; j = 403;
			if (i <= tempmesh.GetNE() && j <= tempmesh.GetNSE())
			if (tempmesh.SurfaceElement(j).GetIndex()==2)
			{
			const Element & el = tempmesh.VolumeElement(i);
			const Element2d & sel = tempmesh.SurfaceElement(j);
		
			const Point3d *tripp[3];
			const Point3d *pp[4];
			int tetpi[4], tripi[3];
		
			for (k = 1; k <= 4; k++)
			{
			pp[k-1] = &tempmesh.Point(el.PNum(k));
			tetpi[k-1] = el.PNum(k);
			}
		
			for (k = 1; k <= 3; k++)
			{
			tripp[k-1] = &tempmesh.Point (sel.PNum(k));
			tripi[k-1] = sel.PNum(k);
			}
		
			(*testout) << "Check Triangle " << j << ":";
			for (k = 1; k <= 3; k++)
			(*testout) << " " << sel.PNum(k);
			for (k = 1; k <= 3; k++)
			(*testout) << " " << tempmesh.Point(sel.PNum(k));
			(*testout) << endl;
		
			(*testout) << "Check Tet " << i << ":";
			for (k = 1; k <= 4; k++)
			(*testout) << " " << el.PNum(k);
			for (k = 1; k <= 4; k++)
			(*testout) << " " << tempmesh.Point(el.PNum(k));
			(*testout) << endl;
		
			if (IntersectTetTriangle (&pp[0], &tripp[0], tetpi, tripi))
			{
			cout << "Intesection detected !!" << endl;
			}
			}
		
			tempmesh.Save ("temp.vol");
		
			// end bug search
			*/
        
        
			for (int i = ne; i >= 1; i--)
			{
			if (outer.Test(i))
			{
			  tempels.DeleteElement(i);
			}
			}
        
        
			// mesh.points.SetSize(mesh.points.Size()-4);
        
			for (int i = 0; i < tempels.Size(); i++)
			{
			Element el = new Element(4);
			for (int j = 0; j < 4; j++)
			{
			  el[j] = tempels[i][j];
			}
			mesh.AddVolumeElement(el);
			}
        
			PrintMessage(5, "outer removed");
        
			mesh.FindOpenElements(domainnr);
        
			mesh.Compress();
			PopStatus();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void LoadRules(string filename, string[] prules)
		{
		  string buf = new string(new char[256]);
		  istream ist;
		  string tr1 = null;
        
		  if (filename)
		  {
			  PrintMessage(3, "rule-filename = ", filename);
			  ist = new ifstream(filename);
		  }
		  else
		  {
			  /* connect tetrules to one string */
			  PrintMessage(3, "Use internal rules");
			  if (!prules)
			  {
				  prules = tetrules;
			  }
        
		//C++ TO C# CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged:
			  string * hcp = prules;
			  uint len = 0;
			  while (*hcp)
			  {
			  len += Convert.ToStringhcp.Length;
			  hcp = hcp.Substring(1);
			  }
			  tr1 = new string(new char[len]);
			  tr1 = null;
			  hcp = prules; //  tetrules;
        
        
		//C++ TO C# CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged:
			  char * tt1 = tr1;
			  while (*hcp)
			  {
			  tt1 += *hcp;
			  tt1 += Convert.ToStringhcp.Length;
			  hcp = hcp.Substring(1);
			  }
        
        
		#if WIN32
			  // VC++ 2005 workaround
			  for (uint i = 0; i < len; i++)
			  {
			if (tr1[i] == ',')
			{
			  tr1[i] = ':';
			}
			  }
		#endif
        
			  ist = new istringstream(tr1);
		  }
        
		  if (!ist.good())
		  {
			  cerr << "Rule description file " << filename << " not found" << "\n";
			  ist = null;
			  Environment.Exit(1);
		  }
        
		  while (!ist.eof())
		  {
			  buf = null;
			  ist >> buf;
        
			  if (string.Compare(buf, "rule") == 0)
			  {
			  vnetrule rule = new vnetrule();
			  rule.LoadRule(ist);
			  rules.Append(rule);
			  if (rule.TestOk() == 0)
			  {
				  PrintSysError("Parser3d: Rule ", rules.Size(), " not ok");
				  Environment.Exit(1);
			  }
			  }
			  else if (string.Compare(buf, "tolfak") == 0)
			  {
			  ist >> tolfak;
			  }
		  }
		  ist = null;
		  tr1 = null;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int ApplyRules(Array<Point3d, PointIndex.BASE> lpoints, Array<int, PointIndex.BASE> allowpoint, Array<MiniElement2d> lfaces, int lfacesplit, INDEX_2_HASHTABLE<int> connectedpairs, Array<Element> elements, Array<int> delfaces, int tolerance, double sloppy, int rotind1, ref float retminerr)
        
		{
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//  static Timer t("ruler3 - all");
		  RegionTimer reg = new RegionTimer(ApplyRules_t);
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//  static Timer tstart("ruler3 - rule start");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//  static Timer tloop("ruler3 - rule loop");
        
		  ApplyRules_tstart.Start();
		  float err;
		  float minerr;
		  float teterr;
		  float minteterr;
		  char ok;
		  char found;
		  char hc;
		  // vnetrule * rule;
		  Vector oldu = new Vector();
		  Vector newu = new Vector();
		  Vector newu1 = new Vector();
		  Vector newu2 = new Vector();
		  Vector allp = new Vector();
		  Vec3d ui = new Vec3d();
		  Point3d np = new Point3d();
		  MiniElement2d locface = null;
		  int loktestmode;
        
        
		  Array<int, PointIndex.BASE> pused = new Array<int, PointIndex.BASE>(); // point is already mapped, number of uses
		  ArrayMem<char,100> fused = new ArrayMem<char,100>(); // face is already mapped
		  ArrayMem<PointIndex,100> pmap = new ArrayMem<PointIndex,100>(); // map of reference point to local point
		  ArrayMem<bool,100> pfixed = new ArrayMem<bool,100>(); // point mapped by face-map
		  ArrayMem<int,100> fmapi = new ArrayMem<int,100>(); // face in reference is mapped to face nr ...
		  ArrayMem<int,100> fmapr = new ArrayMem<int,100>(); // face in reference is rotated to map
		  ArrayMem<Point3d,100> transfreezone = new ArrayMem<Point3d,100>(); // transformed free-zone
		  INDEX_2_CLOSED_HASHTABLE<int> ledges = new INDEX_2_CLOSED_HASHTABLE<int>(100); // edges in local environment
        
		  ArrayMem<Point3d,100> tempnewpoints = new ArrayMem<Point3d,100>();
		  Array<MiniElement2d> tempnewfaces = new Array<MiniElement2d>();
		  ArrayMem<int,100> tempdelfaces = new ArrayMem<int,100>();
		  Array<Element> tempelements = new Array<Element>();
		  ArrayMem<Box3d,100> triboxes = new ArrayMem<Box3d,100>(); // bounding boxes of local faces
        
		  Array<int, PointIndex.BASE> pnearness = new Array<int, PointIndex.BASE>();
		  Array<int> fnearness = new Array<int>();
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//  static int cnt = 0;
		  ApplyRules_cnt++;
        
		  delfaces.SetSize(0);
		  elements.SetSize(0);
        
		  // determine topological distance of faces and points to
		  // base element
        
		  pnearness.SetSize(lpoints.Size());
		  fnearness.SetSize(lfacesplit);
        
		  pnearness = INT_MAX / 10;
        
		  foreach (PointIndex pi in lfaces[0].PNums())
		  {
			pnearness[pi] = 0;
		  }
        
		  NgProfiler.RegionTimer reg2 = new NgProfiler.RegionTimer(98);
        
		  NgProfiler.StartTimer(90);
        
		  for (int loop = 0; loop < 2; loop++)
		  {
        
			  for (int i = 0; i < lfacesplit; i++)
			  {
			  MiniElement2d hface = lfaces[i];
        
			  int minn = INT_MAX - 1;
			  foreach (PointIndex pi in hface.PNums())
			  {
				  int hi = pnearness[pi];
				  if (hi < minn)
				  {
					  minn = hi;
				  }
			  }
			  if (minn < INT_MAX / 10)
			  {
				foreach (PointIndex pi in hface.PNums())
				{
				  if (pnearness[pi] > minn + 1)
				  {
				pnearness[pi] = minn + 1;
				  }
				}
			  }
			  }
        
			  for (int i = 1; i <= connectedpairs.GetNBags(); i++)
			  {
			for (int j = 1; j <= connectedpairs.GetBagSize(i); j++)
			{
				INDEX_2 edge = new INDEX_2();
				int val;
				connectedpairs.GetData(i, j, ref edge, ref val);
        
				if (pnearness[edge.I1()] > pnearness[edge.I2()] + 1)
				{
				  pnearness[edge.I1()] = pnearness[edge.I2()] + 1;
				}
        
				if (pnearness[edge.I2()] > pnearness[edge.I1()] + 1)
				{
				  pnearness[edge.I2()] = pnearness[edge.I1()] + 1;
				}
			}
			  }
		  }
        
		  foreach (int i in fnearness.Range())
		  {
			  int sum = 0;
			  foreach (PointIndex pi in lfaces[i].PNums())
			  {
				sum += pnearness[pi];
			  }
			  fnearness[i] = sum;
		  }
        
        
		  NgProfiler.StopTimer(90);
		  NgProfiler.StartTimer(91);
        
		  // find bounding boxes of faces
        
		  triboxes.SetSize(lfaces.Size());
		  // for (int i = 0; i < lfaces.Size(); i++)
		  foreach (var i in lfaces.Range())
		  {
			  MiniElement2d face = lfaces[i];
			  triboxes[i].SetPoint(lpoints[face[0]]);
			  for (int j = 1; j < face.GetNP(); j++)
			  {
			triboxes[i].AddPoint(lpoints[face[j]]);
			  }
		  }
        
		  NgProfiler.StopTimer(91);
		  NgProfiler.StartTimer(92);
        
        
		  bool useedges = false;
		  for (int ri = 0; ri < rules.Size(); ri++)
		  {
			if (rules[ri].GetNEd())
			{
				useedges = true;
			}
		  }
        
		  if (useedges)
		  {
			  ledges.SetSize((uint)(5 * lfacesplit));
        
			  for (int j = 0; j < lfacesplit; j++)
			  {
			// if (fnearness[j] <= 5) 
				MiniElement2d face = lfaces[j];
				int newp;
				int oldp;
        
				newp = face[face.GetNP() - 1];
				for (int k = 0; k < face.GetNP(); k++)
				{
				oldp = newp;
				newp = face[k];
				ledges.Set(INDEX_2.Sort(oldp, newp), 1);
				}
			  }
		  }
        
		  NgProfiler.StopTimer(92);
        
		  NgProfiler.RegionTimer reg3 = new NgProfiler.RegionTimer(99);
        
		  pused.SetSize(lpoints.Size());
		  fused.SetSize(lfaces.Size());
        
		  found = 0;
		  minerr = tolfak * tolerance * tolerance;
		  minteterr = sloppy * tolerance;
        
		  if (testmode)
		  {
			(*testout) << "cnt = " << ApplyRules_cnt << " class = " << tolerance << "\n";
		  }
        
        
        
		  // impossible, if no rule can be applied at any tolerance class
		  bool impossible = true;
        
        
		  // check each rule:
		  ApplyRules_tstart.Stop();
		  ApplyRules_tloop.Start();
		  for (int ri = 1; ri <= rules.Size(); ri++)
		  {
			  int @base = (lfaces[0].GetNP() == 3) ? 100 : 200;
			  NgProfiler.RegionTimer regx1 = new NgProfiler.RegionTimer(@base);
			  NgProfiler.RegionTimer regx = new NgProfiler.RegionTimer(@base + ri);
        
			  // sprintf (problems.Elem(ri), "");
			  *problems.Elem(ri) = '\0';
        
			  vnetrule rule = rules.Get(ri);
        
			  if (rule.GetNP(1) != lfaces[0].GetNP())
			  {
			continue;
			  }
        
			  if (rule.GetQuality() > tolerance)
			  {
			  if (rule.GetQuality() < 100)
			  {
				  impossible = false;
			  }
        
			  if (testmode)
			  {
				problems.Elem(ri) = "Quality not ok";
			  }
			  continue;
			  }
        
			  if (testmode)
			  {
			problems.Elem(ri) = "no mapping found";
			  }
        
			  loktestmode = testmode || rule.TestFlag('t') || tolerance > 5;
        
			  if (loktestmode != 0)
			  {
			(*testout) << "Rule " << ri << " = " << rule.Name() << "\n";
			  }
        
			  pmap.SetSize(rule.GetNP());
			  fmapi.SetSize(rule.GetNF());
			  fmapr.SetSize(rule.GetNF());
        
			  fused = 0;
			  pused = 0;
			  foreach (var p in pmap)
			  {
				  p.Invalidate();
			  }
			  fmapi = 0;
        
			  foreach (int i in fmapr.Range())
			  {
				fmapr[i] = rule.GetNP(i + 1);
			  }
        
			  fused[0] = 1;
			  fmapi[0] = 1;
			  fmapr[0] = rotind1;
        
			  for (int j = 1; j <= lfaces[0].GetNP(); j++)
			  {
			  PointIndex locpi = lfaces[0].PNumMod(j + rotind1);
			  pmap.Set(rule.GetPointNr(1, j), locpi);
			  pused[locpi]++;
			  }
        
			  /*
			map all faces
			nfok .. first nfok-1 faces are mapped properly
			*/
        
			  int nfok = 2;
			  NgProfiler.RegionTimer regfa = new NgProfiler.RegionTimer(300);
			  NgProfiler.RegionTimer regx2 = new NgProfiler.RegionTimer(@base+50 + ri);
			  while (nfok >= 2)
			  {
        
			  if (nfok <= rule.GetNOldF())
			  {
				  // not all faces mapped
        
				  ok = 0;
				  int locfi = fmapi.Get(nfok);
				  int locfr = fmapr.Get(nfok);
        
				  int actfnp = rule.GetNP(nfok);
        
				  while (!ok)
				  {
				  locfr++;
				  if (locfr == actfnp + 1)
				  {
					  locfr = 1;
					  locfi++;
					  if (locfi > lfacesplit)
					  {
						  break;
					  }
				  }
        
        
				  if (fnearness.Get(locfi) > rule.GetFNearness(nfok) || fused.Get(locfi) || actfnp != lfaces.Get(locfi).GetNP())
				  {
					  // face not feasible in any rotation
        
					  locfr = actfnp;
				  }
				  else
				  {
        
					  ok = 1;
        
					  locface = lfaces.Get(locfi);
        
        
					  // reference point already mapped differently ?
					  for (int j = 1; j <= actfnp && ok; j++)
					  {
					  PointIndex locpi = pmap.Get(rule.GetPointNr(nfok, j));
					  if (locpi.IsValid() && locpi != locface.PNumMod(j + locfr))
					  {
						ok = 0;
					  }
					  }
        
					  // local point already used or point outside tolerance ?
					  for (int j = 1; j <= actfnp && ok; j++)
					  {
					  int refpi = rule.GetPointNr(nfok, j);
        
					  if (!pmap.Get(refpi).IsValid())
					  {
						  PointIndex locpi = locface.PNumMod(j + locfr);
        
						  if (pused[locpi])
						  {
						ok = 0;
						  }
						  else
						  {
						  Point3d lp = lpoints[locpi];
						  Point3d rp = rule.GetPoint(refpi);
        
						  if (Dist2(lp, rp) * rule.PointDistFactor(refpi) > minerr)
						  {
							  impossible = false;
							  ok = 0;
						  }
						  }
					  }
					  }
				  }
				  }
        
        
				  if (ok)
				  {
				  // map face nfok
        
				  fmapi.Set(nfok, locfi);
				  fmapr.Set(nfok, locfr);
				  fused.Set(locfi, 1);
        
				  for (int j = 1; j <= rule.GetNP(nfok); j++)
				  {
					  PointIndex locpi = locface.PNumMod(j + locfr);
        
					  if (rule.GetPointNr(nfok, j) <= 3 && pmap.Get(rule.GetPointNr(nfok, j)) != locpi)
					  {
					(*testout) << "change face1 point, mark1" << "\n";
					  }
        
					  pmap.Set(rule.GetPointNr(nfok, j), locpi);
					  pused[locpi]++;
				  }
        
				  nfok++;
				  }
				  else
				  {
				  // backtrack one face
				  fmapi.Set(nfok, 0);
				  fmapr.Set(nfok, rule.GetNP(nfok));
				  nfok--;
        
				  fused.Set(fmapi.Get(nfok), 0);
				  for (int j = 1; j <= rule.GetNP(nfok); j++)
				  {
					  int refpi = rule.GetPointNr(nfok, j);
					  pused[pmap.Get(refpi)]--;
        
					  if (pused[pmap.Get(refpi)] == 0)
					  {
					  // pmap.Set(refpi, 0);
								  pmap.Elem(refpi).Invalidate();
					  }
				  }
				  }
			  }
        
			  else
        
			  {
				  NgProfiler.RegionTimer regfb = new NgProfiler.RegionTimer(301);
        
				  // all faces are mapped
				  // now map all isolated points:
        
				  if (loktestmode != 0)
				  {
				  (*testout) << "Faces Ok" << "\n";
				  problems.Elem(ri) = "Faces Ok";
				  }
        
				  int npok = 1;
				  int incnpok = 1;
        
				  pfixed.SetSize(pmap.Size());
					  /*
				  for (int i = 1; i <= pmap.Size(); i++)
				pfixed.Set(i, (pmap.Get(i) != 0) );
					  */
					  foreach (int i in pmap.Range())
					  {
						pfixed[i] = pmap[i].IsValid();
					  }
        
				  while (npok >= 1)
				  {
        
				  if (npok <= rule.GetNOldP())
				  {
        
					  if (pfixed.Get(npok))
        
					  {
					  if (incnpok != 0)
					  {
						npok++;
					  }
					  else
					  {
						npok--;
					  }
					  }
        
					  else
        
					  {
					  PointIndex locpi = pmap.Elem(npok);
					  ok = 0;
        
					  if (locpi.IsValid())
					  {
						pused[locpi]--;
					  }
        
					  while (!ok && locpi < lpoints.Size() - 1 + PointIndex.BASE)
					  {
						  ok = 1;
						  locpi++;
        
						  if (pused[locpi] || pnearness[locpi] > rule.GetPNearness(npok))
						  {
						  ok = 0;
						  }
						  else if (allowpoint[locpi] != 2)
						  {
						  ok = 0;
						  if (allowpoint[locpi] == 1)
						  {
							impossible = false;
						  }
						  }
						  else
						  {
						  Point3d lp = lpoints[locpi];
						  Point3d rp = rule.GetPoint(npok);
        
						  if (Dist2(lp, rp) * rule.PointDistFactor(npok) > minerr)
						  {
							  ok = 0;
							  impossible = false;
						  }
						  }
					  }
        
        
					  if (ok)
					  {
						  pmap.Set(npok, locpi);
        
						  if (npok <= 3)
						  {
						(*testout) << "set face1 point, mark3" << "\n";
						  }
        
						  pused[locpi]++;
						  npok++;
						  incnpok = 1;
					  }
        
					  else
        
					  {
						  // pmap.Set (npok, 0);
									  pmap.Elem(npok).Invalidate();
        
						  if (npok <= 3)
						  {
						(*testout) << "set face1 point, mark4" << "\n";
						  }
        
						  npok--;
						  incnpok = 0;
					  }
					  }
				  }
        
				  else
        
				  {
					  NgProfiler.RegionTimer regfa2 = new NgProfiler.RegionTimer(302);
        
					  // all points are mapped
        
					  if (loktestmode != 0)
					  {
					  (*testout) << "Mapping found!!: Rule " << rule.Name() << "\n";
					  foreach (var pi in pmap)
					  {
						(*testout) << pi << " ";
					  }
					  (*testout) << "\n";
					  problems.Elem(ri) = "mapping found";
					  (*testout) << rule.GetNP(1) << " = " << lfaces.Get(1).GetNP() << "\n";
					  }
        
					  ok = 1;
        
        
					  // check mapedges:
					  for (int i = 1; i <= rule.GetNEd(); i++)
					  {
					  INDEX_2 in2 = new INDEX_2(pmap.Get(rule.GetEdge(i).i1), pmap.Get(rule.GetEdge(i).i2));
					  in2.Sort();
					  if (!ledges.Used(in2))
					  {
						  ok = 0;
					  }
					  }
        
        
					  // check prism edges:
					  for (int i = 1; i <= rule.GetNE(); i++)
					  {
					  Element el = rule.GetElement(i);
					  if (el.GetType() == ELEMENT_TYPE.PRISM)
					  {
						  for (int j = 1; j <= 3; j++)
						  {
						  INDEX_2 in2 = new INDEX_2(pmap.Get(el.PNum(j)), pmap.Get(el.PNum(j + 3)));
						  in2.Sort();
						  if (!connectedpairs.Used(in2))
						  {
							  ok = 0;
						  }
						  }
					  }
					  if (el.GetType() == ELEMENT_TYPE.PYRAMID)
					  {
						  if (loktestmode != 0)
						  {
						(*testout) << "map pyramid, rule = " << rule.Name() << "\n";
						  }
						  for (int j = 1; j <= 2; j++)
						  {
						  INDEX_2 in2 = new INDEX_2();
						  if (j == 1)
						  {
							  in2.I1() = pmap.Get(el.PNum(2));
							  in2.I2() = pmap.Get(el.PNum(3));
						  }
						  else
						  {
							  in2.I1() = pmap.Get(el.PNum(1));
							  in2.I2() = pmap.Get(el.PNum(4));
						  }
						  in2.Sort();
						  if (!connectedpairs.Used(in2))
						  {
							  ok = 0;
							  if (loktestmode != 0)
							  {
							(*testout) << "no pair" << "\n";
							  }
						  }
						  }
					  }
        
					  }
        
        
        
					  for (int i = rule.GetNOldF() + 1; i <= rule.GetNF(); i++)
					  {
					fmapi.Set(i, 0);
					  }
        
        
					  if (ok)
					  {
					  foundmap.Elem(ri)++;
					  }
        
        
        
        
					  // deviation of existing points
        
					  oldu.SetSize(3 * rule.GetNOldP());
					  newu.SetSize(3 * (rule.GetNP() - rule.GetNOldP()));
					  allp.SetSize(3 * rule.GetNP());
        
					  for (int i = 1; i <= rule.GetNOldP(); i++)
					  {
					  Point3d lp = lpoints[pmap.Get(i)];
					  Point3d rp = rule.GetPoint(i);
					  oldu(3 * i - 3) = lp.X() - rp.X();
								  oldu(3 * i - 2) = lp.Y() - rp.Y();
					  oldu(3 * i - 1) = lp.Z() - rp.Z();
        
					  allp(3 * i - 3) = lp.X();
								  allp(3 * i - 2) = lp.Y();
								  allp(3 * i - 1) = lp.Z();
					  }
        
					  if (rule.GetNP() > rule.GetNOldP())
					  {
					  newu.SetSize(rule.GetOldUToNewU.functorMethod().Height());
					  rule.GetOldUToNewU.functorMethod().Mult(oldu, newu);
					  }
        
					  //		      int idiff = 3 * (rule->GetNP()-rule->GetNOldP());
					  int idiff = 3 * rule.GetNOldP();
					  for (int i = rule.GetNOldP() + 1; i <= rule.GetNP(); i++)
					  {
					  Point3d rp = rule.GetPoint(i);
					  allp(3 * i - 3) = rp.X() + newu(3 * i - 3 - idiff);
								  allp(3 * i - 2) = rp.Y() + newu(3 * i - 2 - idiff);
								  allp(3 * i - 1) = rp.Z() + newu(3 * i - 1 - idiff);
					  }
        
					  rule.SetFreeZoneTransformation(allp, tolerance + (int)sloppy);
        
					  if (rule.ConvexFreeZone() == 0)
					  {
					  ok = 0;
					  problems.Elem(ri) = "Freezone not convex";
        
					  if (loktestmode != 0)
					  {
						(*testout) << "Freezone not convex" << "\n";
					  }
					  }
        
					  if (loktestmode != 0)
					  {
					  Array<Point3d> fz = rule.GetTransFreeZone();
					  (*testout) << "Freezone: " << "\n";
					  for (int i = 1; i <= fz.Size(); i++)
					  {
						(*testout) << fz.Get(i) << "\n";
					  }
					  }
        
        
					  // check freezone:
        
					  for (int i = 1; i <= lpoints.Size(); i++)
					  {
					  if (!pused.Get(i))
					  {
						  Point3d lp = lpoints.Get(i);
        
						  if (rule.fzbox.IsIn(lp) != 0)
						  {
						  if (rule.IsInFreeZone(lp) != 0)
						  {
							  if (loktestmode != 0)
							  {
							  (*testout) << "Point " << i << " in Freezone" << "\n";
							  problems.Elem(ri) = string.Format("locpoint {0:D} in Freezone", i);
							  }
							  ok = 0;
							  break;
						  }
						  }
					  }
					  }
        
					  for (int i = 1; i <= lfaces.Size() && ok; i++)
					  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//			  static Array<int> lpi(4);
        
					  if (!fused.Get(i))
					  {
						  int triin;
						  MiniElement2d lfacei = lfaces.Get(i);
        
						  if (!triboxes.Elem(i).Intersect(rule.fzbox))
						  {
						triin = 0;
						  }
						  else
						  {
						  int li;
						  int lj;
						  for (li = 1; li <= lfacei.GetNP(); li++)
						  {
							  int lpii = 0;
							  PointIndex pi = lfacei.PNum(li);
							  for (lj = 1; lj <= rule.GetNOldP(); lj++)
							  {
							if (pmap.Get(lj) == pi)
							{
							  lpii = lj;
							}
							  }
							  ApplyRules_lpi.Elem(li) = lpii;
						  }
        
        
						  if (lfacei.GetNP() == 3)
						  {
							  triin = rule.IsTriangleInFreeZone(lpoints[lfacei.PNum(1)], lpoints[lfacei.PNum(2)], lpoints[lfacei.PNum(3)], ApplyRules_lpi, 1);
						  }
						  else
						  {
							  triin = rule.IsQuadInFreeZone(lpoints[lfacei.PNum(1)], lpoints[lfacei.PNum(2)], lpoints[lfacei.PNum(3)], lpoints[lfacei.PNum(4)], ApplyRules_lpi, 1);
						  }
						  }
        
        
						  if (triin == -1)
						  {
						  ok = 0;
						  }
        
						  if (triin == 1)
						  {
		#if TEST_JS
						  ok = 0;
        
						  if (loktestmode != 0)
						  {
							  (*testout) << "El with " << lfaces.Get(i).GetNP() << " points in freezone: " << lfaces.Get(i).PNum(1) << " - " << lfaces.Get(i).PNum(2) << " - " << lfaces.Get(i).PNum(3) << " - " << lfaces.Get(i).PNum(4) << "\n";
							  for (int lj = 1; lj <= lfaces.Get(i).GetNP(); lj++)
							  {
							(*testout) << lpoints[lfaces.Get(i).PNum(lj)] << " ";
							  }
        
							  (*testout) << "\n";
        
							  problems.Elem(ri) = string.Format("triangle ({0:D}, {1:D}, {2:D}) in Freezone", lfaces.Get(i).PNum(1), lfaces.Get(i).PNum(2), lfaces.Get(i).PNum(3));
						  }
		#else
						  if (loktestmode != 0)
						  {
							  if (lfacei.GetNP() == 3)
							  {
							  (*testout) << "Triangle in freezone: " << lfacei.PNum(1) << " - " << lfacei.PNum(2) << " - " << lfacei.PNum(3) << ", or " << lpoints[lfacei.PNum(1)] << " - " << lpoints[lfacei.PNum(2)] << " - " << lpoints[lfacei.PNum(3)] << "\n";
							  (*testout) << "lpi = " << ApplyRules_lpi.Get(1) << ", " << ApplyRules_lpi.Get(2) << ", " << ApplyRules_lpi.Get(3) << "\n";
							  }
							  else
							  {
							  (*testout) << "Quad in freezone: " << lfacei.PNum(1) << " - " << lfacei.PNum(2) << " - " << lfacei.PNum(3) << " - " << lfacei.PNum(4) << ", or " << lpoints[lfacei.PNum(1)] << " - " << lpoints[lfacei.PNum(2)] << " - " << lpoints[lfacei.PNum(3)] << " - " << lpoints[lfacei.PNum(4)] << "\n";
							  }
        
							  problems.Elem(ri) = string.Format("triangle ({0:D}, {1:D}, {2:D}) in Freezone", (int)(lfaces.Get(i).PNum(1)), (int)(lfaces.Get(i).PNum(2)), (int)(lfaces.Get(i).PNum(3)));
						  }
        
						  hc = 0;
						  for (int k = rule.GetNOldF() + 1; k <= rule.GetNF(); k++)
						  {
							  if (rule.GetPointNr(k, 1) <= rule.GetNOldP() && rule.GetPointNr(k, 2) <= rule.GetNOldP() && rule.GetPointNr(k, 3) <= rule.GetNOldP())
							  {
							  for (int j = 1; j <= 3; j++)
							  {
								if (lfaces.Get(i).PNumMod(j) == pmap.Get(rule.GetPointNr(k, 1)) && lfaces.Get(i).PNumMod(j + 1) == pmap.Get(rule.GetPointNr(k, 3)) && lfaces.Get(i).PNumMod(j + 2) == pmap.Get(rule.GetPointNr(k, 2)))
								{
								fmapi.Elem(k) = i;
								hc = 1;
        
        
		 // 						(*testout) << "found from other side: " 
		//  							   << rule->Name() 
		//  							   << " ( " << pmap.Get (rule->GetPointNr(k, 1))
		//  							   << " - " << pmap.Get (rule->GetPointNr(k, 2))
		//  							   << " - " << pmap.Get (rule->GetPointNr(k, 3)) << " ) "
		//  							   << endl;
        
								problems.Elem(ri) = "other";
								}
							  }
							  }
						  }
        
						  if (!hc)
						  {
							  if (loktestmode != 0)
							  {
							  (*testout) << "Triangle in freezone: " << lfaces.Get(i).PNum(1) << " - " << lfaces.Get(i).PNum(2) << " - " << lfaces.Get(i).PNum(3) << "\n";
        
							  problems.Elem(ri) = string.Format("triangle ({0:D}, {1:D}, {2:D}) in Freezone", (int)(lfaces.Get(i).PNum(1)), (int)(lfaces.Get(i).PNum(2)), (int)(lfaces.Get(i).PNum(3)));
							  }
							  ok = 0;
						  }
		#endif
						  }
					  }
        
					  }
        
        
					  if (ok)
					  {
					  err = 0F;
					  for (int i = 1; i <= rule.GetNOldP(); i++)
					  {
						  double hf = rule.CalcPointDist(i, lpoints[pmap.Get(i)]);
						  if (hf > err)
						  {
							  err = hf;
						  }
					  }
        
        
					  if (loktestmode != 0)
					  {
						  (*testout) << "Rule ok" << "\n";
						  problems.Elem(ri) = string.Format("Rule ok, err = {0:f}", err);
					  }
        
        
					  //			  newu = rule->GetOldUToNewU() * oldu;
        
					  // set new points:
								  int oldnp = rule.GetNOldP();
					  int noldlp = lpoints.Size();
					  int noldlf = lfaces.Size();
        
					  for (int i = oldnp + 1; i <= rule.GetNP(); i++)
					  {
						  np.CopyFrom(rule.GetPoint(i));
						  np.X() += newu(3 * (i - oldnp) - 3);
						  np.Y() += newu(3 * (i - oldnp) - 2);
						  np.Z() += newu(3 * (i - oldnp) - 1);
						  lpoints.Append(np);
									  pmap.Elem(i) = lpoints.Size() - 1 + PointIndex.BASE;
					  }
        
					  // Set new Faces:
        
					  for (int i = rule.GetNOldF() + 1; i <= rule.GetNF(); i++)
					  {
						if (!fmapi.Get(i))
						{
						MiniElement2d nface = new MiniElement2d(rule.GetNP(i));
						for (int j = 1; j <= nface.GetNP(); j++)
						{
						  nface.PNum(j) = pmap.Get(rule.GetPointNr(i, j));
						}
        
						lfaces.Append(nface);
						}
					  }
        
					  // Delete old Faces:
        
					  for (int i = 1; i <= rule.GetNDelF(); i++)
					  {
						delfaces.Append(fmapi.Get(rule.GetDelFace(i)));
					  }
					  for (int i = rule.GetNOldF() + 1; i <= rule.GetNF(); i++)
					  {
						if (fmapi.Get(i))
						{
						delfaces.Append(fmapi.Get(i));
						fmapi.Elem(i) = 0;
						}
					  }
        
        
					  // check orientation
					  for (int i = 1; i <= rule.GetNO() && ok; i++)
					  {
						  fourint fouri;
        
						  fouri = rule.GetOrientation(i);
						  Vec3d v1 = new Vec3d(lpoints[pmap.Get(fouri.i1)], lpoints[pmap.Get(fouri.i2)]);
						  Vec3d v2 = new Vec3d(lpoints[pmap.Get(fouri.i1)], lpoints[pmap.Get(fouri.i3)]);
						  Vec3d v3 = new Vec3d(lpoints[pmap.Get(fouri.i1)], lpoints[pmap.Get(fouri.i4)]);
        
						  Vec3d n = new Vec3d();
						  Cross(v1, v2, n);
						  //if (n * v3 >= -1e-7*n.Length()*v3.Length()) // OR -1e-7???
						  if (n * v3 >= -1e-9 != null)
						  {
						  if (loktestmode != 0)
						  {
							  problems.Elem(ri) = "Orientation wrong";
							  (*testout) << "Orientation wrong (" << n * v3 << ")" << "\n";
						  }
						  ok = 0;
						  }
					  }
        
        
        
					  // new points in free-zone ?
					  for (int i = rule.GetNOldP() + 1; i <= rule.GetNP() && ok; i++)
					  {
						if (rule.IsInFreeZone(lpoints.Get(pmap.Get(i))) == 0)
						{
						if (loktestmode != 0)
						{
							(*testout) << "Newpoint " << lpoints.Get(pmap.Get(i)) << " outside convex hull" << "\n";
							problems.Elem(ri) = "newpoint outside convex hull";
						}
						ok = 0;
        
						}
					  }
        
					  // insert new elements
        
					  for (int i = 1; i <= rule.GetNE(); i++)
					  {
						  elements.Append(rule.GetElement(i));
						  for (int j = 1; j <= elements.Get(i).NP(); j++)
						  {
						elements.Elem(i).PNum(j) = pmap.Get(elements.Get(i).PNum(j));
						  }
					  }
        
        
					  // Calculate Element badness
        
					  teterr = 0F;
					  for (int i = 1; i <= elements.Size(); i++)
					  {
						  double hf = CalcElementBadness(lpoints, elements.Get(i));
						  if (hf > teterr)
						  {
							  teterr = hf;
						  }
					  }
        
					  /*
					    // keine gute Erfahrung am 25.1.2000, js
					  if (ok && teterr < 100 &&
					      (rule->TestFlag('b') || tolerance > 10) )
					    {
					      (*mycout) << "Reset teterr " 
						   << rule->Name() 
						   << " err = " << teterr 
						   << endl;
					      teterr = 1;
					    }
					  */
        
					  // compare edgelength
					  if (rule.TestFlag('l') != 0)
					  {
						  double oldlen = 0;
						  double newlen = 0;
        
						  for (int i = 1; i <= rule.GetNDelF(); i++)
						  {
						  Element2d face = rule.GetFace(rule.GetDelFace(i));
						  for (int j = 1; j <= 3; j++)
						  {
							  Point3d p1 = lpoints[pmap.Get(face.PNumMod(j))];
							  Point3d p2 = lpoints[pmap.Get(face.PNumMod(j + 1))];
							  oldlen += Dist(p1, p2);
						  }
						  }
        
						  for (int i = rule.GetNOldF() + 1; i <= rule.GetNF(); i++)
						  {
						  Element2d face = rule.GetFace(i);
						  for (int j = 1; j <= 3; j++)
						  {
							  Point3d p1 = lpoints[pmap.Get(face.PNumMod(j))];
							  Point3d p2 = lpoints[pmap.Get(face.PNumMod(j + 1))];
							  newlen += Dist(p1, p2);
						  }
						  }
        
						  if (oldlen < newlen)
						  {
						  ok = 0;
						  if (loktestmode != 0)
						  {
							problems.Elem(ri) = "oldlen < newlen";
						  }
						  }
					  }
        
        
					  if (loktestmode != 0)
					  {
						(*testout) << "ok = " << (int)ok << "teterr = " << teterr << "minteterr = " << minteterr << "\n";
					  }
        
        
					  if (ok && teterr < tolerance)
					  {
						  canuse.Elem(ri)++;
						  /*
						  (*testout) << "can use rule " << rule->Name() 
							 << ", err = " << teterr << endl;
						  for (i = 1; i <= pmap.Size(); i++)
						(*testout) << pmap.Get(i) << " ";
						  (*testout) << endl;
						  */
        
						  if (string.Compare(problems.Elem(ri), "other") == 0)
						  {
						  if (teterr < minother)
						  {
							minother = teterr;
						  }
						  }
						  else
						  {
						  if (teterr < minwithoutother)
						  {
							minwithoutother = teterr;
						  }
						  }
					  }
        
        
					  if (teterr > minteterr)
					  {
						  impossible = false;
					  }
        
					  if (ok && teterr < minteterr)
					  {
        
						  if (loktestmode != 0)
						  {
						(*testout) << "use rule" << "\n";
						  }
        
						  found = ri;
						  minteterr = teterr;
        
						  if (testmode)
						  {
						  for (int i = 1; i <= rule.GetNOldP(); i++)
						  {
							  (*testout) << "P" << i << ": Ref: " << rule.GetPoint(i) << "  is: " << lpoints.Get(pmap.Get(i)) << "\n";
						  }
						  }
        
						  tempnewpoints.SetSize(0);
						  for (int i = noldlp + 1; i <= lpoints.Size(); i++)
						  {
						tempnewpoints.Append(lpoints.Get(i));
						  }
        
						  tempnewfaces.SetSize(0);
						  for (int i = noldlf + 1; i <= lfaces.Size(); i++)
						  {
						tempnewfaces.Append(lfaces.Get(i));
						  }
        
						  tempdelfaces.SetSize(0);
						  for (int i = 1; i <= delfaces.Size(); i++)
						  {
						tempdelfaces.Append(delfaces.Get(i));
						  }
        
						  tempelements.SetSize(0);
						  for (int i = 1; i <= elements.Size(); i++)
						  {
						tempelements.Append(elements.Get(i));
						  }
					  }
        
        
					  lpoints.SetSize(noldlp);
					  lfaces.SetSize(noldlf);
					  delfaces.SetSize(0);
					  elements.SetSize(0);
					  }
        
					  npok = rule.GetNOldP();
					  incnpok = 0;
				  }
				  }
        
				  nfok = rule.GetNOldF();
        
				  for (int j = 1; j <= rule.GetNP(nfok); j++)
				  {
				  int refpi = rule.GetPointNr(nfok, j);
				  pused[pmap.Get(refpi)]--;
        
				  if (pused[pmap.Get(refpi)] == 0)
				  {
							pmap.Elem(refpi).Invalidate();
				  }
				  }
        
			  }
			  }
			  if (loktestmode != 0)
			  {
			(*testout) << "end rule" << "\n";
			  }
		  }
		  ApplyRules_tloop.Stop();
        
		  if (found)
		  {
			  /*
			  for (i = 1; i <= tempnewpoints.Size(); i++)
			lpoints.Append (tempnewpoints.Get(i));
			  */
			  foreach (Point3d p in tempnewpoints)
			  {
				lpoints.Append(p);
			  }
			  /*
			  for (i = 1; i <= tempnewfaces.Size(); i++)
			if (tempnewfaces.Get(i).PNum(1))
			  lfaces.Append (tempnewfaces.Get(i));
			  */
			  foreach (int i in tempnewfaces.Range())
			  {
			if (tempnewfaces[i].PNum(1).IsValid())
			{
			  lfaces.Append(tempnewfaces[i]);
			}
			  }
			  /*
			  for (i = 1; i <= tempdelfaces.Size(); i++)
			delfaces.Append (tempdelfaces.Get(i));
			  */
			  foreach (int i in tempdelfaces.Range())
			  {
				delfaces.Append(tempdelfaces[i]);
			  }
			  /*
			  for (i = 1; i <= tempelements.Size(); i++)
			elements.Append (tempelements.Get(i));
			  */
			  foreach (int i in tempelements.Range())
			  {
				elements.Append(tempelements[i]);
			  }
		  }
        
		  retminerr = minerr;
        
        
		  if (impossible && found == 0)
		  {
			return -1;
		  }
        
		  return found;
		}
	}
}