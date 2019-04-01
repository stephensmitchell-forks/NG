using System;

///
public class MeshOptimize3d
{
  private readonly MeshingParameters mp;
  public MeshOptimize3d(MeshingParameters amp)
  {
	  this.mp = amp;
	  ;
  }

  /*
    Combine two points to one.
    Set new point into the center, if both are
    inner points.
    Connect inner point to boundary point, if one
    point is inner point.
  */


//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer CombineImprove_t("MeshOptimize3d::CombineImprove");

  public void CombineImprove(Mesh mesh, OPTIMIZEGOAL goal = OPT_QUALITY)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("MeshOptimize3d::CombineImprove");
	RegionTimer reg = new RegionTimer(CombineImprove_t);

	int np = mesh.GetNP();
	int ne = mesh.GetNE();

	TABLE<ElementIndex, PointIndex.BASE> elementsonnode = new TABLE<ElementIndex, PointIndex.BASE>(np);
	Array<ElementIndex> hasonepi = new Array<ElementIndex>();
	Array<ElementIndex> hasbothpi = new Array<ElementIndex>();

	Array<double> oneperr = new Array<double>();
	Array<double> elerrs = new Array<double>(ne);

	PrintMessage(3, "CombineImprove");
	(*testout) << "Start CombineImprove" << "\n";

	//  mesh.CalcSurfacesOfNode ();
	string savetask = multithread.task;
	multithread.task = "Combine Improve";


	double totalbad = 0;
	for (ElementIndex ei = 0; ei < ne; ei++)
	{
		if (mesh.GetDimension() == 3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
		{
	  continue;
		}
		double elerr = CalcBad(mesh.Points(), mesh[ei], 0);
		totalbad += elerr;
		elerrs[ei] = elerr;
	}

	if (goal == OPT_QUALITY)
	{
		totalbad = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
		(*testout) << "Total badness = " << totalbad << "\n";
		PrintMessage(5, "Total badness = ", totalbad);
	}

	for (ElementIndex ei = 0; ei < ne; ei++)
	{
	  if (!mesh[ei].IsDeleted())
	  {
		for (int j = 0; j < mesh[ei].GetNP(); j++)
		{
	  elementsonnode.Add(mesh[ei][j], ei);
		}
	  }
	}

	INDEX_2_HASHTABLE<int> edgetested = new INDEX_2_HASHTABLE<int>(np + 1);

	int cnt = 0;

	for (ElementIndex ei = 0; ei < ne; ei++)
	{
		if (mesh.GetDimension() == 3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
		{
			continue;
		}
		if (multithread.terminate)
		{
	  break;
		}

		multithread.percent = 100.0 * (ei + 1) / ne;

		if (mesh.ElementType(ei) == FIXEDELEMENT)
		{
	  continue;
		}

		for (int j = 0; j < 6; j++)
		{
		Element elemi = mesh[ei];
		if (elemi.IsDeleted())
		{
			continue;
		}

		int[][] tetedges =
		{
			new int[] {0, 1},
			new int[] {0, 2},
			new int[] {0, 3},
			new int[] {1, 2},
			new int[] {1, 3},
			new int[] {2, 3}
		};

		PointIndex pi1 = elemi[tetedges[j][0]];
		PointIndex pi2 = elemi[tetedges[j][1]];

		if (pi2 < pi1)
		{
			netgen.GlobalMembers.Swap(ref pi1, ref pi2);
		}

		INDEX_2 si2 = new INDEX_2(pi1, pi2);
		si2.Sort();

		if (edgetested.Used(si2))
		{
			continue;
		}
		edgetested.Set(si2, 1);


		// hasonepoint.SetSize(0);
		//	  hasbothpoints.SetSize(0);
		hasonepi.SetSize(0);
		hasbothpi.SetSize(0);

		FlatArray<ElementIndex> row1 = elementsonnode[pi1];
		for (int k = 0; k < row1.Size(); k++)
		{
			Element elem = mesh[row1[k]];
			if (elem.IsDeleted())
			{
				continue;
			}

			if (elem[0] == pi2 || elem[1] == pi2 || elem[2] == pi2 || elem[3] == pi2)
			{
			hasbothpi.Append(row1[k]);
			}
			else
			{
			hasonepi.Append(row1[k]);
			}
		}

		FlatArray<ElementIndex> row2 = elementsonnode[pi2];
		for (int k = 0; k < row2.Size(); k++)
		{
			Element elem = mesh[row2[k]];
			if (elem.IsDeleted())
			{
				continue;
			}

			if (elem[0] == pi1 || elem[1] == pi1 || elem[2] == pi1 || elem[3] == pi1)
			{
		  ;
			}
			else
			{
			hasonepi.Append(row2[k]);
			}
		}

		double bad1 = 0;
		for (int k = 0; k < hasonepi.Size(); k++)
		{
		  bad1 += elerrs[hasonepi[k]];
		}
		for (int k = 0; k < hasbothpi.Size(); k++)
		{
		  bad1 += elerrs[hasbothpi[k]];
		}

		MeshPoint p1 = mesh[pi1];
		MeshPoint p2 = mesh[pi2];


		// if (mesh.PointType(pi2) != INNERPOINT)
		if (p2.Type() != INNERPOINT)
		{
		  continue;
		}

		MeshPoint pnew = new MeshPoint();
		// if (mesh.PointType(pi1) != INNERPOINT)
		if (p1.Type() != INNERPOINT)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pnew = p1;
		  pnew.CopyFrom(p1);
		}
		else
		{
		  pnew = netgen.GlobalMembers.Center(p1, p2);
		}

		mesh[pi1] = pnew;
		mesh[pi2] = pnew;

		oneperr.SetSize(hasonepi.Size());

		double bad2 = 0;
		for (int k = 0; k < hasonepi.Size(); k++)
		{
			Element elem = mesh[hasonepi[k]];
			double err = CalcBad(mesh.Points(), elem, 0);
			// CalcTetBadness (mesh[elem[0]], mesh[elem[1]],
			// mesh[elem[2]], mesh[elem[3]], 0, mparam);
			bad2 += err;
			oneperr[k] = err;
		}

		mesh[pi1] = p1;
		mesh[pi2] = p2;

		// if (mesh.PointType(pi1) != INNERPOINT)
		if (p1.Type() != INNERPOINT)
		{
			for (int k = 0; k < hasonepi.Size(); k++)
			{
			Element elem = mesh[hasonepi[k]];
			int l;
			for (l = 0; l < 4; l++)
			{
			  if (elem[l] == pi2)
			  {
			  elem[l] = pi1;
			  break;
			  }
			}

			elem.flags.illegal_valid = 0;
			if (!mesh.LegalTet(elem))
			{
			  bad2 += 1e4;
			}

			if (l < 4)
			{
				elem.flags.illegal_valid = 0;
				elem[l] = pi2;
			}
			}
		}

		if (bad2 / hasonepi.Size() < bad1 / (hasonepi.Size() + hasbothpi.Size()))
		{
			mesh[pi1] = pnew;
			cnt++;

			FlatArray<ElementIndex> row = elementsonnode[pi2];
			for (int k = 0; k < row.Size(); k++)
			{
			Element elem = mesh[row[k]];
			if (elem.IsDeleted())
			{
				continue;
			}

			elementsonnode.Add(pi1, row[k]);
			for (int l = 0; l < elem.GetNP(); l++)
			{
			  if (elem[l] == pi2)
			  {
				elem[l] = pi1;
			  }
			}

			elem.flags.illegal_valid = 0;
			if (!mesh.LegalTet(elem))
			{
			  (*testout) << "illegal tet " << elementsonnode[pi2][k] << "\n";
			}
			}

			for (int k = 0; k < hasonepi.Size(); k++)
			{
		  elerrs[hasonepi[k]] = oneperr[k];
			}

			for (int k = 0; k < hasbothpi.Size(); k++)
			{
			mesh[hasbothpi[k]].flags.illegal_valid = 0;
			mesh[hasbothpi[k]].Delete();
			}
		}
		}
	}

	mesh.Compress();
	mesh.MarkIllegalElements();

	PrintMessage(5, cnt, " elements combined");
	(*testout) << "CombineImprove done" << "\n";

	totalbad = 0;
	for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	{
	  if (!(mesh.GetDimension() == 3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex()))
	  {
		totalbad += CalcBad(mesh.Points(), mesh[ei], 0);
	  }
	}

	if (goal == OPT_QUALITY)
	{
		totalbad = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
		(*testout) << "Total badness = " << totalbad << "\n";

		int cntill = 0;
		for (ElementIndex ei = 0; ei < ne; ei++)
		{
	  if (!(mesh.GetDimension() == 3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex()))
	  {
		if (!mesh.LegalTet(mesh[ei]))
		{
		  cntill++;
		}
	  }
		}

		PrintMessage(5, cntill, " illegal tets");
	}
	multithread.task = savetask;
  }


  /*
    Mesh improvement by edge splitting.
    If mesh quality is improved by inserting a node into an inner edge,
    the edge is split into two parts.
  */
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer SplitImprove_t("MeshOptimize3d::SplitImprove");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer SplitImprove_tloop("MeshOptimize3d::SplitImprove loop");

  public void SplitImprove(Mesh mesh, OPTIMIZEGOAL goal = OPT_QUALITY)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("MeshOptimize3d::SplitImprove");
	RegionTimer reg = new RegionTimer(SplitImprove_t);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer tloop("MeshOptimize3d::SplitImprove loop");

	double bad1;
	double bad2;
	double badmax;
	double badlimit;
	int cnt = 0;
	int np = mesh.GetNP();
	int ne = mesh.GetNE();

	TABLE<ElementIndex,PointIndex.BASE> elementsonnode = new TABLE<ElementIndex,PointIndex.BASE>(np);
	Array<ElementIndex> hasbothpoints = new Array<ElementIndex>();

	BitArray origpoint = new BitArray(np + 1); // big enough for 0 and 1-based
	BitArray boundp = new BitArray(np + 1);
	origpoint.Set();

	Array<double> elerrs = new Array<double>(ne);
	BitArray illegaltet = new BitArray(ne);
	illegaltet.Clear();

	string savetask = multithread.task;
	multithread.task = "Split Improve";

	PrintMessage(3, "SplitImprove");
	(*testout) << "start SplitImprove" << "\n";

	Array<INDEX_3> locfaces = new Array<INDEX_3>();

	INDEX_2_HASHTABLE<int> edgetested = new INDEX_2_HASHTABLE<int>(np);

	bad1 = 0;
	badmax = 0;
	for (ElementIndex ei = 0; ei < ne; ei++)
	{
		if (mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
		{
			continue;
		}

		elerrs[ei] = CalcBad(mesh.Points(), mesh[ei], 0);
		bad1 += elerrs[ei];
		if (elerrs[ei] > badmax)
		{
			badmax = elerrs[ei];
		}
	}

	PrintMessage(5, "badmax = ", badmax);
	badlimit = 0.5 * badmax;

	boundp.Clear();
	foreach (var el in mesh.SurfaceElements())
	{
	  foreach (PointIndex pi in el.PNums())
	  {
		boundp.Set(pi);
	  }
	}

	if (goal == OPT_QUALITY)
	{
		bad1 = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
		(*testout) << "Total badness = " << bad1 << "\n";
	}

	foreach (ElementIndex ei in mesh.VolumeElements().Range())
	{
	  foreach (PointIndex pi in mesh[ei].PNums())
	  {
		elementsonnode.Add(pi, ei);
	  }
	}

	mesh.MarkIllegalElements();
	if (goal == OPT_QUALITY || goal == OPT_LEGAL)
	{
		int cntill = 0;
		for (ElementIndex ei = 0; ei < ne; ei++)
		{
			//	  if (!LegalTet (volelements.Get(i)))
			if (mesh[ei].flags.illegal)
			{
				cntill++;
				illegaltet.Set(ei);
			}
		}
	}

	SplitImprove_tloop.Start();
	foreach (ElementIndex ei in mesh.VolumeElements().Range())
	{
		Element elem = mesh[ei];

		if (mp.only3D_domain_nr && mp.only3D_domain_nr != elem.GetIndex())
		{
			continue;
		}
		if (multithread.terminate)
		{
	  break;
		}

		multithread.percent = 100.0 * (ei + 1) / ne;

		bool ltestmode = false;

		if (elerrs[ei] < badlimit && !illegaltet.Test(ei))
		{
			continue;
		}

		if ((goal == OPT_LEGAL) && !illegaltet.Test(ei) && CalcBad(mesh.Points(), elem, 0) < 1e3)
		{
	  continue;
		}

		if (ltestmode)
		{
		(*testout) << "test el " << ei << "\n";
		for (int j = 0; j < 4; j++)
		{
		  (*testout) << elem[j] << " ";
		}
		(*testout) << "\n";
		}

		for (int j = 0; j < 6; j++)
		{
		int[][] tetedges =
		{
			new int[] {0, 1},
			new int[] {0, 2},
			new int[] {0, 3},
			new int[] {1, 2},
			new int[] {1, 3},
			new int[] {2, 3}
		};

		PointIndex pi1 = elem[tetedges[j][0]];
		PointIndex pi2 = elem[tetedges[j][1]];

		if (pi2 < pi1)
		{
			netgen.GlobalMembers.Swap(ref pi1, ref pi2);
		}
		if (pi2 >= elementsonnode.Size() + PointIndex.BASE)
		{
			continue; // old number of points
		}

		if (!origpoint.Test(pi1) || !origpoint.Test(pi2))
		{
		  continue;
		}

		INDEX_2 i2 = new INDEX_2(pi1, pi2);
		i2.Sort();

		if (mesh.BoundaryEdge(pi1, pi2))
		{
			continue;
		}

		if (edgetested.Used(i2) && !illegaltet.Test(ei))
		{
			continue;
		}
		edgetested.Set(i2, 1);

		hasbothpoints.SetSize(0);
			/*
		for (int k = 1; k <= elementsonnode.EntrySize(pi1); k++)
		  {
			ElementIndex elnr = elementsonnode.Get(pi1, k);
			*/
			foreach (ElementIndex ei in elementsonnode[pi1])
			{
			Element el = mesh[ei];
				bool has1 = el.PNums().Contains(pi1);
				bool has2 = el.PNums().Contains(pi2);

			if (has1 && has2)
			{
		  if (!hasbothpoints.Contains(ei))
		  {
			hasbothpoints.Append(ei);
		  }
			}
			}

		bad1 = 0;

			foreach (ElementIndex ei in hasbothpoints)
			{
		  bad1 += CalcBad(mesh.Points(), mesh[ei], 0);
			}

		bool puretet = true;
		foreach (ElementIndex ei in hasbothpoints)
		{
		  if (mesh[ei].GetType() != TET)
		  {
			puretet = false;
		  }
		}
		if (!puretet)
		{
			continue;
		}

		Point3d p1 = mesh[pi1];
		Point3d p2 = mesh[pi2];

		/*
		  pnew = Center (p1, p2);
		
		  points.Elem(pi1) = pnew;
		  bad2 = 0;
		  for (k = 1; k <= hasbothpoints.Size(); k++)
		  bad2 += CalcBad (points,
		  volelements.Get(hasbothpoints.Get(k)), 0);
  
		  points.Elem(pi1) = p1;
		  points.Elem(pi2) = pnew;
		
		  for (k = 1; k <= hasbothpoints.Size(); k++)
		  bad2 += CalcBad (points,
		  volelements.Get(hasbothpoints.Get(k)), 0);
		  points.Elem(pi2) = p2;
		*/

		locfaces.SetSize(0);
		foreach (ElementIndex ei in hasbothpoints)
		{
			Element el = mesh[ei];

			for (int l = 0; l < 4; l++)
			{
		  if (el[l] == pi1 || el[l] == pi2)
		  {
			  INDEX_3 i3 = new INDEX_3();
			  Element2d face = new Element2d(TRIG);
			  el.GetFace(l + 1, face);
			  for (int kk = 1; kk <= 3; kk++)
			  {
				i3.I(kk) = face.PNum(kk);
			  }
			  locfaces.Append(i3);
		  }
			}
		}

		PointFunction1 pf = new PointFunction1(mesh.Points(), locfaces, mp, -1);
		OptiParameters par = new OptiParameters();
		par.maxit_linsearch = 50;
		par.maxit_bfgs = 20;

		Point3d pnew = netgen.GlobalMembers.Center(p1, p2);
		Vector px = new Vector(3);
		px(0) = pnew.X();
		px(1) = pnew.Y();
		px(2) = pnew.Z();

		if (elerrs[ei] > 0.1 * badmax)
		{
		  BFGS(px, pf, par);
		}

		bad2 = pf.Func(px);

		pnew.X() = px(0);
		pnew.Y() = px(1);
		pnew.Z() = px(2);

		PointIndex hpinew = mesh.AddPoint(pnew);
		//	  ptyps.Append (INNERPOINT);

		for (int k = 0; k < hasbothpoints.Size(); k++)
		{
			Element oldel = mesh[hasbothpoints[k]];
			Element newel1 = new Element(oldel);
			Element newel2 = new Element(oldel);

			oldel.flags.illegal_valid = 0;
			newel1.flags.illegal_valid = 0;
			newel2.flags.illegal_valid = 0;

			for (int l = 0; l < 4; l++)
			{
			if (newel1[l] == pi2)
			{
				newel1[l] = hpinew;
			}
			if (newel2[l] == pi1)
			{
				newel2[l] = hpinew;
			}
			}

			if (!mesh.LegalTet(oldel))
			{
				bad1 += 1e6;
			}
			if (!mesh.LegalTet(newel1))
			{
				bad2 += 1e6;
			}
			if (!mesh.LegalTet(newel2))
			{
				bad2 += 1e6;
			}
		}

		// mesh.PointTypes().DeleteLast();
		mesh.Points().DeleteLast();

		if (bad2 < bad1)
		  /*	      (bad1 > 1e4 && boundp.Test(pi1) && boundp.Test(pi2)) */
		{
			cnt++;

			PointIndex pinew = mesh.AddPoint(pnew);

			foreach (ElementIndex ei in hasbothpoints)
			{
			Element oldel = mesh[ei];
			Element newel = new Element(oldel);

			newel.flags.illegal_valid = 0;
			oldel.flags.illegal_valid = 0;

			for (int l = 0; l < 4; l++)
			{
				origpoint.Clear(oldel[l]);

				if (oldel[l] == pi2)
				{
					oldel[l] = pinew;
				}
				if (newel[l] == pi1)
				{
					newel[l] = pinew;
				}
			}
			mesh.AddVolumeElement(newel);
			}

			j = 10; // end j-loop
		}
		}
	}
	SplitImprove_tloop.Stop();

	mesh.Compress();
	PrintMessage(5, cnt, " splits performed");

	(*testout) << "Splitt - Improve done" << "\n";

	if (goal == OPT_QUALITY)
	{
		bad1 = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
		(*testout) << "Total badness = " << bad1 << "\n";

		int cntill = 0;
		ne = mesh.GetNE();
		for (ElementIndex ei = 0; ei < ne; ei++)
		{
		  if (!mesh.LegalTet(mesh[ei]))
		  {
			cntill++;
		  }
		}
		//      cout << cntill << " illegal tets" << endl;
	}

	multithread.task = savetask;
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer SwapImprove_t("MeshOptimize3d::SwapImprove");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer SwapImprove_tloop("MeshOptimize3d::SwapImprove loop");

  public void SwapImprove(Mesh mesh, OPTIMIZEGOAL goal = OPT_QUALITY, BitArray working_elements = null)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("MeshOptimize3d::SwapImprove");
	RegionTimer reg = new RegionTimer(SwapImprove_t);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer tloop("MeshOptimize3d::SwapImprove loop");

	PointIndex pi3 = new PointIndex(0);
	PointIndex pi4 = new PointIndex(0);
	PointIndex pi5 = new PointIndex(0);
	PointIndex pi6 = new PointIndex(0);
	int cnt = 0;

	Element el21 = new Element(TET);
	Element el22 = new Element(TET);
	Element el31 = new Element(TET);
	Element el32 = new Element(TET);
	Element el33 = new Element(TET);
	Element el1 = new Element(TET);
	Element el2 = new Element(TET);
	Element el3 = new Element(TET);
	Element el4 = new Element(TET);
	Element el1b = new Element(TET);
	Element el2b = new Element(TET);
	Element el3b = new Element(TET);
	Element el4b = new Element(TET);

	double bad1;
	double bad2;
	double bad3;

	int np = mesh.GetNP();
	int ne = mesh.GetNE();

	// contains at least all elements at node
	TABLE<ElementIndex,PointIndex.BASE> elementsonnode = new TABLE<ElementIndex,PointIndex.BASE>(np);

	Array<ElementIndex> hasbothpoints = new Array<ElementIndex>();

	PrintMessage(3, "SwapImprove ");
	(*testout) << "\n" << "Start SwapImprove" << "\n";

	string savetask = multithread.task;
	multithread.task = "Swap Improve";

	//  mesh.CalcSurfacesOfNode ();

	INDEX_3_HASHTABLE<int> faces = new INDEX_3_HASHTABLE<int>(mesh.GetNOpenElements() / 3 + 2);
	if (goal == OPT_CONFORM)
	{
		for (int i = 1; i <= mesh.GetNOpenElements(); i++)
		{
		Element2d hel = mesh.OpenElement(i);
		INDEX_3 face = new INDEX_3(hel[0], hel[1], hel[2]);
		face.Sort();
		faces.Set(face, 1);
		}
	}

	// Calculate total badness
	if (goal == OPT_QUALITY)
	{
		bad1 = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
		(*testout) << "Total badness = " << bad1 << "\n";
	}

	// find elements on node
	for (ElementIndex ei = 0; ei < ne; ei++)
	{
	  foreach (PointIndex pi in mesh[ei].PNums())
	  {
		elementsonnode.Add(pi, ei);
	  }
	}
	  /*
	  for (int j = 0; j < mesh[ei].GetNP(); j++)
	    elementsonnode.Add (mesh[ei][j], ei);
	  */


	// INDEX_2_HASHTABLE<int> edgeused(2 * ne + 5);
	INDEX_2_CLOSED_HASHTABLE<int> edgeused = new INDEX_2_CLOSED_HASHTABLE<int>((uint)(12 * ne + 5));

	SwapImprove_tloop.Start();
	for (ElementIndex ei = 0; ei < ne; ei++)
	{
		if (multithread.terminate)
		{
	  break;
		}

		if (mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
		{
			continue;
		}

		multithread.percent = 100.0 * (ei + 1) / ne;

		if ((mesh.ElementType(ei)) == FIXEDELEMENT)
		{
	  continue;
		}

		if (working_elements != null && ei < working_elements.Size() && !working_elements.Test(ei))
		{
	  continue;
		}

		if (mesh[ei].IsDeleted())
		{
	  continue;
		}

		if ((goal == OPT_LEGAL) && mesh.LegalTet(mesh[ei]) && CalcBad(mesh.Points(), mesh[ei], 0) < 1e3)
		{
	  continue;
		}

		//      int onlybedges = 1;

		for (int j = 0; j < 6; j++)
		{
		// loop over edges

		Element elemi = mesh[ei];
		if (elemi.IsDeleted())
		{
			continue;
		}

		int mattyp = elemi.GetIndex();

		int[][] tetedges =
		{
			new int[] {0, 1},
			new int[] {0, 2},
			new int[] {0, 3},
			new int[] {1, 2},
			new int[] {1, 3},
			new int[] {2, 3}
		};

		PointIndex pi1 = elemi[tetedges[j][0]];
		PointIndex pi2 = elemi[tetedges[j][1]];

		if (pi2 < pi1)
		{
			netgen.GlobalMembers.Swap(ref pi1, ref pi2);
		}

		if (mesh.BoundaryEdge(pi1, pi2))
		{
			continue;
		}


		INDEX_2 i2 = new INDEX_2(pi1, pi2);
		i2.Sort();
		if (edgeused.Used(i2))
		{
			continue;
		}
		edgeused.Set(i2, 1);

		hasbothpoints.SetSize(0);
		for (int k = 0; k < elementsonnode[pi1].Size(); k++)
		{
			bool has1 = false;
			bool has2 = false;
			ElementIndex elnr = elementsonnode[pi1][k];
			Element elem = mesh[elnr];

			if (elem.IsDeleted())
			{
				continue;
			}

			for (int l = 0; l < elem.GetNP(); l++)
			{
			if (elem[l] == pi1)
			{
				has1 = true;
			}
			if (elem[l] == pi2)
			{
				has2 = true;
			}
			}

			if (has1 && has2)
			{ // only once
					if (hasbothpoints.Contains(elnr))
					{
					  has1 = false;
					}

			if (has1)
			{
						hasbothpoints.Append(elnr);
			}
			}
		}

		bool puretet = true;
		foreach (ElementIndex ei in hasbothpoints)
		{
		  if (mesh[ei].GetType() != TET)
		  {
			puretet = false;
		  }
		}

		if (!puretet)
		{
			continue;
		}

		int nsuround = hasbothpoints.Size();

		if (nsuround == 3)
		{
			Element elem = mesh[hasbothpoints[0]];
			for (int l = 0; l < 4; l++)
			{
		  if (elem[l] != pi1 && elem[l] != pi2)
		  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi4 = pi3;
			  pi4.CopyFrom(pi3);
			  pi3 = elem[l];
		  }
			}

			el31[0] = pi1;
			el31[1] = pi2;
			el31[2] = pi3;
			el31[3] = pi4;
			el31.SetIndex(mattyp);

			if (WrongOrientation(mesh.Points(), el31))
			{
			netgen.GlobalMembers.Swap(ref pi3, ref pi4);
			el31[2] = pi3;
			el31[3] = pi4;
			}

			pi5 = 0;
			for (int k = 0; k < 3; k++) // JS, 201212
			{
			Element elemk = mesh[hasbothpoints[k]];
			bool has1 = false;
			for (int l = 0; l < 4; l++)
			{
			  if (elemk[l] == pi4)
			  {
				has1 = true;
			  }
			}
			if (has1)
			{
				for (int l = 0; l < 4; l++)
				{
			  if (elemk[l] != pi1 && elemk[l] != pi2 && elemk[l] != pi4)
			  {
				pi5 = elemk[l];
			  }
				}
			}
			}

			if (!pi5.IsValid())
			{
				  throw new Exception("Illegal state observed in SwapImprove");
			}


			el32[0] = pi1;
			el32[1] = pi2;
			el32[2] = pi4;
			el32[3] = pi5;
			el32.SetIndex(mattyp);

			el33[0] = pi1;
			el33[1] = pi2;
			el33[2] = pi5;
			el33[3] = pi3;
			el33.SetIndex(mattyp);

			elementsonnode.Add(pi4, hasbothpoints[1]);
			elementsonnode.Add(pi3, hasbothpoints[2]);

			bad1 = CalcBad(mesh.Points(), el31, 0) + CalcBad(mesh.Points(), el32, 0) + CalcBad(mesh.Points(), el33, 0);

			el31.flags.illegal_valid = 0;
			el32.flags.illegal_valid = 0;
			el33.flags.illegal_valid = 0;

			if (!mesh.LegalTet(el31) || !mesh.LegalTet(el32) || !mesh.LegalTet(el33))
			{
		  bad1 += 1e4;
			}

			el21[0] = pi3;
			el21[1] = pi4;
			el21[2] = pi5;
			el21[3] = pi2;
			el21.SetIndex(mattyp);

			el22[0] = pi5;
			el22[1] = pi4;
			el22[2] = pi3;
			el22[3] = pi1;
			el22.SetIndex(mattyp);

			bad2 = CalcBad(mesh.Points(), el21, 0) + CalcBad(mesh.Points(), el22, 0);

			el21.flags.illegal_valid = 0;
			el22.flags.illegal_valid = 0;

			if (!mesh.LegalTet(el21) || !mesh.LegalTet(el22))
			{
		  bad2 += 1e4;
			}


			if (goal == OPT_CONFORM && bad2 < 1e4)
			{
			INDEX_3 face = new INDEX_3(pi3, pi4, pi5);
			face.Sort();
			if (faces.Used(face))
			{
				// (*testout) << "3->2 swap, could improve conformity, bad1 = " << bad1
				//				 << ", bad2 = " << bad2 << endl;
				if (bad2 < 1e4)
				{
			  bad1 = 2 * bad2;
				}
			}
			/*
			  else
			  {
			  INDEX_2 hi1(pi3, pi4);
			  hi1.Sort();
			  INDEX_2 hi2(pi3, pi5);
			  hi2.Sort();
			  INDEX_2 hi3(pi4, pi5);
			  hi3.Sort();
  
			  if (boundaryedges->Used (hi1) ||
			  boundaryedges->Used (hi2) ||
			  boundaryedges->Used (hi3) )
			  bad1 = 2 * bad2;
			  }
			*/
			}

			if (bad2 < bad1)
			{
			//		  (*mycout) << "3->2 " << flush;
			//		  (*testout) << "3->2 conversion" << endl;
			cnt++;


			/*
			(*testout) << "3->2 swap, old els = " << endl
				   << mesh[hasbothpoints[0]] << endl
				   << mesh[hasbothpoints[1]] << endl
				   << mesh[hasbothpoints[2]] << endl
				   << "new els = " << endl
				   << el21 << endl
				   << el22 << endl;
			*/

			el21.flags.illegal_valid = 0;
			el22.flags.illegal_valid = 0;
			mesh[hasbothpoints[0]] = el21;
			mesh[hasbothpoints[1]] = el22;
			for (int l = 0; l < 4; l++)
			{
			  mesh[hasbothpoints[2]][l].Invalidate();
			}
			mesh[hasbothpoints[2]].Delete();

			for (int k = 0; k < 2; k++)
			{
			  for (int l = 0; l < 4; l++)
			  {
				elementsonnode.Add(mesh[hasbothpoints[k]][l], hasbothpoints[k]);
			  }
			}
			}
		}

		if (nsuround == 4)
		{
			Element elem1 = mesh[hasbothpoints[0]];
			for (int l = 0; l < 4; l++)
			{
		  if (elem1[l] != pi1 && elem1[l] != pi2)
		  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi4 = pi3;
			  pi4.CopyFrom(pi3);
			  pi3 = elem1[l];
		  }
			}

			el1[0] = pi1;
			el1[1] = pi2;
			el1[2] = pi3;
			el1[3] = pi4;
			el1.SetIndex(mattyp);

			if (WrongOrientation(mesh.Points(), el1))
			{
			netgen.GlobalMembers.Swap(ref pi3, ref pi4);
			el1[2] = pi3;
			el1[3] = pi4;
			}

			pi5.Invalidate();
			for (int k = 0; k < 4; k++)
			{
			Element elem = mesh[hasbothpoints[k]];
					bool has1 = elem.PNums().Contains(pi4);
			if (has1)
			{
				for (int l = 0; l < 4; l++)
				{
			  if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi4)
			  {
				pi5 = elem[l];
			  }
				}
			}
			}

			pi6.Invalidate();
			for (int k = 0; k < 4; k++)
			{
			Element elem = mesh[hasbothpoints[k]];
					bool has1 = elem.PNums().Contains(pi3);
			if (has1)
			{
				for (int l = 0; l < 4; l++)
				{
			  if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi3)
			  {
				pi6 = elem[l];
			  }
				}
			}
			}


			/*
			INDEX_2 i22(pi3, pi5);
			i22.Sort();
			INDEX_2 i23(pi4, pi6);
			i23.Sort();
			*/

			el1[0] = pi1;
			el1[1] = pi2;
			el1[2] = pi3;
			el1[3] = pi4;
			el1.SetIndex(mattyp);

			el2[0] = pi1;
			el2[1] = pi2;
			el2[2] = pi4;
			el2[3] = pi5;
			el2.SetIndex(mattyp);

			el3[0] = pi1;
			el3[1] = pi2;
			el3[2] = pi5;
			el3[3] = pi6;
			el3.SetIndex(mattyp);

			el4[0] = pi1;
			el4[1] = pi2;
			el4[2] = pi6;
			el4[3] = pi3;
			el4.SetIndex(mattyp);

			//        elementsonnode.Add (pi4, hasbothpoints.Elem(2));
			//        elementsonnode.Add (pi3, hasbothpoints.Elem(3));

			bad1 = CalcBad(mesh.Points(), el1, 0) + CalcBad(mesh.Points(), el2, 0) + CalcBad(mesh.Points(), el3, 0) + CalcBad(mesh.Points(), el4, 0);


			el1.flags.illegal_valid = 0;
			el2.flags.illegal_valid = 0;
			el3.flags.illegal_valid = 0;
			el4.flags.illegal_valid = 0;


			if (goal != OPT_CONFORM)
			{
			if (!mesh.LegalTet(el1) || !mesh.LegalTet(el2) || !mesh.LegalTet(el3) || !mesh.LegalTet(el4))
			{
			  bad1 += 1e4;
			}
			}

			el1[0] = pi3;
			el1[1] = pi5;
			el1[2] = pi2;
			el1[3] = pi4;
			el1.SetIndex(mattyp);

			el2[0] = pi3;
			el2[1] = pi5;
			el2[2] = pi4;
			el2[3] = pi1;
			el2.SetIndex(mattyp);

			el3[0] = pi3;
			el3[1] = pi5;
			el3[2] = pi1;
			el3[3] = pi6;
			el3.SetIndex(mattyp);

			el4[0] = pi3;
			el4[1] = pi5;
			el4[2] = pi6;
			el4[3] = pi2;
			el4.SetIndex(mattyp);

			bad2 = CalcBad(mesh.Points(), el1, 0) + CalcBad(mesh.Points(), el2, 0) + CalcBad(mesh.Points(), el3, 0) + CalcBad(mesh.Points(), el4, 0);

			el1.flags.illegal_valid = 0;
			el2.flags.illegal_valid = 0;
			el3.flags.illegal_valid = 0;
			el4.flags.illegal_valid = 0;

			if (goal != OPT_CONFORM)
			{
			if (!mesh.LegalTet(el1) || !mesh.LegalTet(el2) || !mesh.LegalTet(el3) || !mesh.LegalTet(el4))
			{
			  bad2 += 1e4;
			}
			}


			el1b[0] = pi4;
			el1b[1] = pi6;
			el1b[2] = pi3;
			el1b[3] = pi2;
			el1b.SetIndex(mattyp);

			el2b[0] = pi4;
			el2b[1] = pi6;
			el2b[2] = pi2;
			el2b[3] = pi5;
			el2b.SetIndex(mattyp);

			el3b[0] = pi4;
			el3b[1] = pi6;
			el3b[2] = pi5;
			el3b[3] = pi1;
			el3b.SetIndex(mattyp);

			el4b[0] = pi4;
			el4b[1] = pi6;
			el4b[2] = pi1;
			el4b[3] = pi3;
			el4b.SetIndex(mattyp);

			bad3 = CalcBad(mesh.Points(), el1b, 0) + CalcBad(mesh.Points(), el2b, 0) + CalcBad(mesh.Points(), el3b, 0) + CalcBad(mesh.Points(), el4b, 0);

			el1b.flags.illegal_valid = 0;
			el2b.flags.illegal_valid = 0;
			el3b.flags.illegal_valid = 0;
			el4b.flags.illegal_valid = 0;

			if (goal != OPT_CONFORM)
			{
			if (!mesh.LegalTet(el1b) || !mesh.LegalTet(el2b) || !mesh.LegalTet(el3b) || !mesh.LegalTet(el4b))
			{
			  bad3 += 1e4;
			}
			}


			/*
			int swap2 = (bad2 < bad1) && (bad2 < bad3);
			int swap3 = !swap2 && (bad3 < bad1);
			
			if ( ((bad2 < 10 * bad1) ||
			  (bad2 < 1e6)) && mesh.BoundaryEdge (pi3, pi5))
		  swap2 = 1;
			else  if ( ((bad3 < 10 * bad1) ||
				(bad3 < 1e6)) && mesh.BoundaryEdge (pi4, pi6))
		  {
			swap3 = 1;
			swap2 = 0;
		  }
			*/
			bool swap2;
			bool swap3;

			if (goal != OPT_CONFORM)
			{
			swap2 = (bad2 < bad1) && (bad2 < bad3);
			swap3 = !swap2 && (bad3 < bad1);
			}
			else
			{
			if (mesh.BoundaryEdge(pi3, pi5))
			{
				bad2 /= 1e6;
			}
			if (mesh.BoundaryEdge(pi4, pi6))
			{
				bad3 /= 1e6;
			}

			swap2 = (bad2 < bad1) && (bad2 < bad3);
			swap3 = !swap2 && (bad3 < bad1);
			}


			if (swap2 || swap3)
			{
			// (*mycout) << "4->4 " << flush;
			cnt++;
			//		  (*testout) << "4->4 conversion" << "\n";
			/*
			  (*testout) << "bad1 = " << bad1
			  << " bad2 = " << bad2
			  << " bad3 = " << bad3 << "\n";
			
			  (*testout) << "Points: " << pi1 << " " << pi2 << " " << pi3
			  << " " << pi4 << " " << pi5 << " " << pi6 << "\n";
			  (*testout) << "Elements: "
			  << hasbothpoints.Get(1) << "  "
			  << hasbothpoints.Get(2) << "  "
			  << hasbothpoints.Get(3) << "  "
			  << hasbothpoints.Get(4) << "  " << "\n";
			*/

			/*
			  {
			  int i1, j1;
			  for (i1 = 1; i1 <= 4; i1++)
			  {
			  for (j1 = 1; j1 <= 4; j1++)
			  (*testout) << volelements.Get(hasbothpoints.Get(i1)).PNum(j1)
			  << "  ";
			  (*testout) << "\n";
			  }
			  }
			*/
			}


			if (swap2)
			{
			//		  (*mycout) << "bad1 = " << bad1 << " bad2 = " << bad2 << "\n";


			/*
			(*testout) << "4->4 swap A, old els = " << endl
				   << mesh[hasbothpoints[0]] << endl
				   << mesh[hasbothpoints[1]] << endl
				   << mesh[hasbothpoints[2]] << endl
				   << mesh[hasbothpoints[3]] << endl
				   << "new els = " << endl
				   << el1 << endl
				   << el2 << endl
				   << el3 << endl
				   << el4 << endl;
			*/



			el1.flags.illegal_valid = 0;
			el2.flags.illegal_valid = 0;
			el3.flags.illegal_valid = 0;
			el4.flags.illegal_valid = 0;

			mesh[hasbothpoints[0]] = el1;
			mesh[hasbothpoints[1]] = el2;
			mesh[hasbothpoints[2]] = el3;
			mesh[hasbothpoints[3]] = el4;

			for (int k = 0; k < 4; k++)
			{
			  for (int l = 0; l < 4; l++)
			  {
				elementsonnode.Add(mesh[hasbothpoints[k]][l], hasbothpoints[k]);
			  }
			}
			}
			else if (swap3)
			{
			// (*mycout) << "bad1 = " << bad1 << " bad3 = " << bad3 << "\n";
			el1b.flags.illegal_valid = 0;
			el2b.flags.illegal_valid = 0;
			el3b.flags.illegal_valid = 0;
			el4b.flags.illegal_valid = 0;


			/*
			(*testout) << "4->4 swap A, old els = " << endl
				   << mesh[hasbothpoints[0]] << endl
				   << mesh[hasbothpoints[1]] << endl
				   << mesh[hasbothpoints[2]] << endl
				   << mesh[hasbothpoints[3]] << endl
				   << "new els = " << endl
				   << el1b << endl
				   << el2b << endl
				   << el3b << endl
				   << el4b << endl;
			*/


			mesh[hasbothpoints[0]] = el1b;
			mesh[hasbothpoints[1]] = el2b;
			mesh[hasbothpoints[2]] = el3b;
			mesh[hasbothpoints[3]] = el4b;

			for (int k = 0; k < 4; k++)
			{
			  for (int l = 0; l < 4; l++)
			  {
				elementsonnode.Add(mesh[hasbothpoints[k]][l], hasbothpoints[k]);
			  }
			}
			}
		}

		if (nsuround >= 5)
		{
			Element hel = new Element(TET);

			ArrayMem<PointIndex, 50> suroundpts = new ArrayMem<PointIndex, 50>((uint)nsuround);
			ArrayMem<bool, 50> tetused = new ArrayMem<bool, 50>((uint)nsuround);

			Element elem = mesh[hasbothpoints[0]];

			for (int l = 0; l < 4; l++)
			{
		  if (elem[l] != pi1 && elem[l] != pi2)
		  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi4 = pi3;
			  pi4.CopyFrom(pi3);
			  pi3 = elem[l];
		  }
			}

			hel[0] = pi1;
			hel[1] = pi2;
			hel[2] = pi3;
			hel[3] = pi4;
			hel.SetIndex(mattyp);

			if (WrongOrientation(mesh.Points(), hel))
			{
			netgen.GlobalMembers.Swap(ref pi3, ref pi4);
			hel[2] = pi3;
			hel[3] = pi4;
			}


			// suroundpts.SetSize (nsuround);
				suroundpts = -17;
			suroundpts[0] = pi3;
			suroundpts[1] = pi4;

			tetused = false;
			tetused[0] = true;

			for (int l = 2; l < nsuround; l++)
			{
			PointIndex oldpi = suroundpts[l - 1];
			PointIndex newpi = new PointIndex();
					newpi.Invalidate();

			for (int k = 0; k < nsuround && !newpi.IsValid(); k++)
			{
			  if (!tetused[k])
			  {
			  Element nel = mesh[hasbothpoints[k]];
			  for (int k2 = 0; k2 < 4 && !newpi.IsValid(); k2++)
			  {
				if (nel[k2] == oldpi)
				{
					newpi = nel[0] + nel[1] + nel[2] + nel[3] - pi1 - pi2 - oldpi;

					tetused[k] = true;
					suroundpts[l] = newpi;
				}
			  }
			  }
			}
			}


			bad1 = 0;
			for (int k = 0; k < nsuround; k++)
			{
			hel[0] = pi1;
			hel[1] = pi2;
			hel[2] = suroundpts[k];
			hel[3] = suroundpts[(k + 1) % nsuround];
			hel.SetIndex(mattyp);

			bad1 += CalcBad(mesh.Points(), hel, 0);
			}

			//  (*testout) << "nsuround = " << nsuround << " bad1 = " << bad1 << endl;


			int bestl = -1;
			int confface = -1;
			int confedge = -1;
			double badopt = bad1;

			for (int l = 0; l < nsuround; l++)
			{
			bad2 = 0;

			for (int k = l + 1; k <= nsuround + l - 2; k++)
			{
				hel[0] = suroundpts[l];
				hel[1] = suroundpts[k % nsuround];
				hel[2] = suroundpts[(k + 1) % nsuround];
				hel[3] = pi2;

				bad2 += CalcBad(mesh.Points(), hel, 0);
				hel.flags.illegal_valid = 0;
				if (!mesh.LegalTet(hel))
				{
					bad2 += 1e4;
				}

				hel[2] = suroundpts[k % nsuround];
				hel[1] = suroundpts[(k + 1) % nsuround];
				hel[3] = pi1;

				bad2 += CalcBad(mesh.Points(), hel, 0);

				hel.flags.illegal_valid = 0;
				if (!mesh.LegalTet(hel))
				{
					bad2 += 1e4;
				}
			}
			// (*testout) << "bad2," << l << " = " << bad2 << endl;

			if (bad2 < badopt)
			{
				bestl = l;
				badopt = bad2;
			}


			if (goal == OPT_CONFORM)
			{
				 // (bad2 <= 100 * bad1 || bad2 <= 1e6))
				bool nottoobad = (bad2 <= bad1) || (bad2 <= 100 * bad1 && bad2 <= 1e18) || (bad2 <= 1e8);

				for (int k = l + 1; k <= nsuround + l - 2; k++)
				{
				INDEX_3 hi3 = new INDEX_3(suroundpts[l], suroundpts[k % nsuround], suroundpts[(k + 1) % nsuround]);
				hi3.Sort();
				if (faces.Used(hi3))
				{
					// (*testout) << "could improve face conformity, bad1 = " << bad1
					// << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
					if (nottoobad)
					{
				  confface = l;
					}
				}
				}

				for (int k = l + 2; k <= nsuround + l - 2; k++)
				{
				if (mesh.BoundaryEdge(suroundpts[l], suroundpts[k % nsuround]))
				{
					/*
					*testout << "could improve edge conformity, bad1 = " << bad1
					 << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
					*/
					if (nottoobad)
					{
				  confedge = l;
					}
				}
				}
			}
			}

			if (confedge != -1)
			{
		  bestl = confedge;
			}
			if (confface != -1)
			{
		  bestl = confface;
			}

			if (bestl != -1)
			{
			// (*mycout) << nsuround << "->" << 2 * (nsuround-2) << " " << flush;
			cnt++;

			for (int k = bestl + 1; k <= nsuround + bestl - 2; k++)
			{
				int k1;

				hel[0] = suroundpts[bestl];
				hel[1] = suroundpts[k % nsuround];
				hel[2] = suroundpts[(k + 1) % nsuround];
				hel[3] = pi2;
				hel.flags.illegal_valid = 0;

				/*
				(*testout) << nsuround << "-swap, new el,top = "
				   << hel << endl;
				*/
				mesh.AddVolumeElement(hel);

				for (k1 = 0; k1 < 4; k1++)
				{
			  elementsonnode.Add(hel[k1], mesh.GetNE() - 1);
				}


				hel[2] = suroundpts[k % nsuround];
				hel[1] = suroundpts[(k + 1) % nsuround];
				hel[3] = pi1;

				/*
				(*testout) << nsuround << "-swap, new el,bot = "
				   << hel << endl;
				*/

				mesh.AddVolumeElement(hel);

				for (k1 = 0; k1 < 4; k1++)
				{
			  elementsonnode.Add(hel[k1], mesh.GetNE() - 1);
				}
			}

			for (int k = 0; k < nsuround; k++)
			{
				Element rel = mesh[hasbothpoints[k]];
				/*
				(*testout) << nsuround << "-swap, old el = "
				   << rel << endl;
				*/
				rel.Delete();
				for (int k1 = 0; k1 < 4; k1++)
				{
			  rel[k1].Invalidate();
				}
			}
			}
		}
		}

		/*
	  if (onlybedges)
	  {
	  (*testout) << "bad tet: "
	  << volelements.Get(i)[0]
	  << volelements.Get(i)[1]
	  << volelements.Get(i)[2]
	  << volelements.Get(i)[3] << "\n";
  
	  if (!mesh.LegalTet (volelements.Get(i)))
	  cerr << "Illegal tet" << "\n";
	  }
		*/
	}
	//  (*mycout) << endl;
	SwapImprove_tloop.Stop();
	/*
	    cout << "edgeused: ";
	    edgeused.PrintMemInfo(cout);
	*/
	PrintMessage(5, cnt, " swaps performed");





	mesh.Compress();

	/*
	if (goal == OPT_QUALITY)
	  {
	    bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
	    //      (*testout) << "Total badness = " << bad1 << endl;
	  }
	*/

	/*
	  for (i = 1; i <= GetNE(); i++)
	  if (volelements.Get(i)[0])
	  if (!mesh.LegalTet (volelements.Get(i)))
	  {
	  cout << "detected illegal tet, 2" << endl;
	  (*testout) << "detected illegal tet1: " << i << endl;
	  }
	*/

	multithread.task = savetask;
  }

  public void SwapImproveSurface(Mesh mesh, OPTIMIZEGOAL goal = OPT_QUALITY, BitArray working_elements = null, Array< Array<int,PointIndex.BASE> > idmaps = null)
  {
	Array< Array<int,PointIndex.BASE> > locidmaps = new Array< Array<int,PointIndex.BASE> >();
	Array< Array<int,PointIndex.BASE> > used_idmaps;

	if (idmaps != null)
	{
	  used_idmaps = idmaps;
	}
	else
	{
		used_idmaps = locidmaps;

		for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
		{
		if (mesh.GetIdentifications().GetType(i) == Identifications.PERIODIC)
		{
			locidmaps.Append(new Array<int,PointIndex.BASE>());
			mesh.GetIdentifications().GetMap(i,*locidmaps.Last(),true);
		}
		}
	}


	PointIndex pi1 = new PointIndex();
	PointIndex pi2 = new PointIndex();
	PointIndex pi3 = new PointIndex();
	PointIndex pi4 = new PointIndex();
	PointIndex pi5 = new PointIndex();
	PointIndex pi6 = new PointIndex();
	PointIndex pi1other = new PointIndex();
	PointIndex pi2other = new PointIndex();
	int cnt = 0;

	//double bad1, bad2, bad3, sbad;
	double bad1;
	double sbad;
	double h;

	int np = mesh.GetNP();
	int ne = mesh.GetNE();
	int nse = mesh.GetNSE();

	int mattype;
	int othermattype;


	// contains at least all elements at node
	TABLE<ElementIndex,PointIndex.BASE> elementsonnode = new TABLE<ElementIndex,PointIndex.BASE>(np);
	TABLE<SurfaceElementIndex,PointIndex.BASE> surfaceelementsonnode = new TABLE<SurfaceElementIndex,PointIndex.BASE>(np);
	TABLE<int,PointIndex.BASE> surfaceindicesonnode = new TABLE<int,PointIndex.BASE>(np);

	Array<ElementIndex> hasbothpoints = new Array<ElementIndex>();
	Array<ElementIndex> hasbothpointsother = new Array<ElementIndex>();

	PrintMessage(3, "SwapImproveSurface ");
	(*testout) << "\n" << "Start SwapImproveSurface" << "\n";

	string savetask = multithread.task;
	multithread.task = "Swap Improve Surface";



	// find elements on node
	for (ElementIndex ei = 0; ei < ne; ei++)
	{
	  for (int j = 0; j < mesh[ei].GetNP(); j++)
	  {
		elementsonnode.Add(mesh[ei][j], ei);
	  }
	}

	for (SurfaceElementIndex sei = 0; sei < nse; sei++)
	{
	  for (int j = 0; j < mesh[sei].GetNP(); j++)
	  {
	  surfaceelementsonnode.Add(mesh[sei][j], sei);
	  if (!surfaceindicesonnode[mesh[sei][j]].Contains(mesh[sei].GetIndex()))
	  {
		surfaceindicesonnode.Add(mesh[sei][j],mesh[sei].GetIndex());
	  }
	  }
	}

	bool periodic;
	int idnum = -1;

	// INDEX_2_HASHTABLE<int> edgeused(2 * ne + 5);
	INDEX_2_CLOSED_HASHTABLE<int> edgeused = new INDEX_2_CLOSED_HASHTABLE<int>((uint)(12 * ne + 5));

	for (ElementIndex ei = 0; ei < ne; ei++)
	{
		if (multithread.terminate)
		{
	  break;
		}

		multithread.percent = 100.0 * (ei + 1) / ne;

		if (mesh.ElementType(ei) == FIXEDELEMENT)
		{
	  continue;
		}

		if (working_elements != null && ei < working_elements.Size() && !working_elements.Test(ei))
		{
	  continue;
		}

		if (mesh[ei].IsDeleted())
		{
	  continue;
		}

		if ((goal == OPT_LEGAL) && mesh.LegalTet(mesh[ei]) && CalcBad(mesh.Points(), mesh[ei], 0) < 1e3)
		{
	  continue;
		}

		Element elemi = mesh[ei];
		//Element elemi = mesh[ei];
		if (elemi.IsDeleted())
		{
			continue;
		}


		mattype = elemi.GetIndex();

		bool swapped = false;

		for (int j = 0; !swapped && j < 6; j++)
		{
		// loop over edges


		int[][] tetedges =
		{
			new int[] {0, 1},
			new int[] {0, 2},
			new int[] {0, 3},
			new int[] {1, 2},
			new int[] {1, 3},
			new int[] {2, 3}
		};

		pi1 = elemi[tetedges[j][0]];
		pi2 = elemi[tetedges[j][1]];


		if (pi2 < pi1)
		{
		  netgen.GlobalMembers.Swap(ref pi1, ref pi2);
		}


		bool found = false;
		for (int k = 0; !found && k < used_idmaps.Size(); k++)
		{
			if (pi2 < used_idmaps[k].Size() + PointIndex.BASE)
			{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi1other = (*(*used_idmaps)[k])[pi1];
			pi1other.CopyFrom((* used_idmaps[k])[pi1]);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi2other = (*(*used_idmaps)[k])[pi2];
			pi2other.CopyFrom((* used_idmaps[k])[pi2]);
			found = (pi1other != 0 && pi2other != 0 && pi1other != pi1 && pi2other != pi2);
			if (found)
			{
			  idnum = k;
			}
			}
		}
		if (found)
		{
		  periodic = true;
		}
		else
		{
			periodic = false;
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi1other = pi1;
			pi1other.CopyFrom(pi1);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pi2other = pi2;
			pi2other.CopyFrom(pi2);
		}



		if (!mesh.BoundaryEdge(pi1, pi2) || mesh.IsSegment(pi1, pi2))
		{
			continue;
		}

		othermattype = -1;


		INDEX_2 i2 = new INDEX_2(pi1, pi2);
		i2.Sort();
		if (edgeused.Used(i2))
		{
			continue;
		}
		edgeused.Set(i2, 1);
		if (periodic)
		{
			i2.I1() = pi1other;
			i2.I2() = pi2other;
			i2.Sort();
			edgeused.Set(i2,1);
		}


		hasbothpoints.SetSize(0);
		hasbothpointsother.SetSize(0);
		for (int k = 0; k < elementsonnode[pi1].Size(); k++)
		{
			bool has1 = false;
			bool has2 = false;
			ElementIndex elnr = elementsonnode[pi1][k];
			Element elem = mesh[elnr];

			if (elem.IsDeleted())
			{
				continue;
			}

			for (int l = 0; l < elem.GetNP(); l++)
			{
			if (elem[l] == pi1)
			{
				has1 = true;
			}
			if (elem[l] == pi2)
			{
				has2 = true;
			}
			}

			if (has1 && has2)
			{
			if (othermattype == -1 && elem.GetIndex() != mattype)
			{
			  othermattype = elem.GetIndex();
			}

			if (elem.GetIndex() == mattype)
			{
				// only once
				for (int l = 0; l < hasbothpoints.Size(); l++)
				{
			  if (hasbothpoints[l] == elnr)
			  {
				has1 = false;
			  }
				}

				if (has1)
				{
			  hasbothpoints.Append(elnr);
				}
			}
			else if (elem.GetIndex() == othermattype)
			{
				// only once
				for (int l = 0; l < hasbothpointsother.Size(); l++)
				{
			  if (hasbothpointsother[l] == elnr)
			  {
				has1 = false;
			  }
				}

				if (has1)
				{
			  hasbothpointsother.Append(elnr);
				}
			}
			else
			{
				Console.Write("problem with domain indices");
				Console.Write("\n");
				(*testout) << "problem: mattype = " << mattype << ", othermattype = " << othermattype << " elem " << elem << " mt " << elem.GetIndex() << "\n" << " pi1 " << pi1 << " pi2 " << pi2 << "\n";
				(*testout) << "hasbothpoints:" << "\n";
				for (int ii = 0; ii < hasbothpoints.Size(); ii++)
				{
			  (*testout) << mesh[hasbothpoints[ii]] << "\n";
				}
				(*testout) << "hasbothpointsother:" << "\n";
				for (int ii = 0; ii < hasbothpointsother.Size(); ii++)
				{
			  (*testout) << mesh[hasbothpointsother[ii]] << "\n";
				}
			}
			}
		}

		if (hasbothpointsother.Size() > 0 && periodic)
		{
		  throw new Exception("SwapImproveSurface: Assumption about interface/periodicity wrong!");
		}

		if (periodic)
		{
			for (int k = 0; k < elementsonnode[pi1other].Size(); k++)
			{
			bool has1 = false;
			bool has2 = false;
			ElementIndex elnr = elementsonnode[pi1other][k];
			Element elem = mesh[elnr];

			if (elem.IsDeleted())
			{
				continue;
			}

			for (int l = 0; l < elem.GetNP(); l++)
			{
				if (elem[l] == pi1other)
				{
					has1 = true;
				}
				if (elem[l] == pi2other)
				{
					has2 = true;
				}
			}

			if (has1 && has2)
			{
				if (othermattype == -1)
				{
			  othermattype = elem.GetIndex();
				}

				// only once
				for (int l = 0; l < hasbothpointsother.Size(); l++)
				{
			  if (hasbothpointsother[l] == elnr)
			  {
				has1 = false;
			  }
				}

				if (has1)
				{
			  hasbothpointsother.Append(elnr);
				}
			}
			}
		}


		//for(k=0; k<hasbothpoints.Size(); k++)
		//  (*testout) << "hasbothpoints["<<k<<"]: " << mesh[hasbothpoints[k]] << endl;


		SurfaceElementIndex sel1 = -1;
		SurfaceElementIndex sel2 = -1;
		SurfaceElementIndex sel1other = -1;
		SurfaceElementIndex sel2other = -1;
		for (int k = 0; k < surfaceelementsonnode[pi1].Size(); k++)
		{
			bool has1 = false;
			bool has2 = false;
			SurfaceElementIndex elnr = surfaceelementsonnode[pi1][k];
			Element2d elem = mesh[elnr];

			if (elem.IsDeleted())
			{
				continue;
			}

			for (int l = 0; l < elem.GetNP(); l++)
			{
			if (elem[l] == pi1)
			{
				has1 = true;
			}
			if (elem[l] == pi2)
			{
				has2 = true;
			}
			}

			if (has1 && has2 && elnr != sel2)
			{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel1 = sel2;
			sel1.CopyFrom(sel2);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel2 = elnr;
			sel2.CopyFrom(elnr);
			}
		}

		if (periodic)
		{
			for (int k = 0; k < surfaceelementsonnode[pi1other].Size(); k++)
			{
			bool has1 = false;
			bool has2 = false;
			SurfaceElementIndex elnr = surfaceelementsonnode[pi1other][k];
			Element2d elem = mesh[elnr];

			if (elem.IsDeleted())
			{
				continue;
			}

			for (int l = 0; l < elem.GetNP(); l++)
			{
				if (elem[l] == pi1other)
				{
					has1 = true;
				}
				if (elem[l] == pi2other)
				{
					has2 = true;
				}
			}

			if (has1 && has2 && elnr != sel2other)
			{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel1other = sel2other;
				sel1other.CopyFrom(sel2other);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel2other = elnr;
				sel2other.CopyFrom(elnr);
			}
			}
		}
		else
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel1other = sel1;
			sel1other.CopyFrom(sel1);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel2other = sel2;
			sel2other.CopyFrom(sel2);
		}

		//(*testout) << "sel1 " << sel1 << " sel2 " << sel2 << " el " << mesh[sel1] << " resp. " << mesh[sel2] << endl;

		PointIndex sp1 = new PointIndex(0);
		PointIndex sp2 = new PointIndex(0);
		PointIndex sp1other = new PointIndex();
		PointIndex sp2other = new PointIndex();
		for (int l = 0; l < mesh[sel1].GetNP(); l++)
		{
		  if (mesh[sel1][l] != pi1 && mesh[sel1][l] != pi2)
		  {
			sp1 = mesh[sel1][l];
		  }
		}
		for (int l = 0; l < mesh[sel2].GetNP(); l++)
		{
		  if (mesh[sel2][l] != pi1 && mesh[sel2][l] != pi2)
		  {
			sp2 = mesh[sel2][l];
		  }
		}

		if (periodic)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sp1other = (*(*used_idmaps)[idnum])[sp1];
			sp1other.CopyFrom((* used_idmaps[idnum])[sp1]);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sp2other = (*(*used_idmaps)[idnum])[sp2];
			sp2other.CopyFrom((* used_idmaps[idnum])[sp2]);

			bool change = false;
			for (int l = 0; !change && l < mesh[sel1other].GetNP(); l++)
			{
		  change = (sp2other == mesh[sel1other][l]);
			}

			if (change)
			{
			SurfaceElementIndex aux = new SurfaceElementIndex(sel1other);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel1other = sel2other;
			sel1other.CopyFrom(sel2other);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sel2other = aux;
			sel2other.CopyFrom(aux);
			}

		}
		else
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sp1other = sp1;
			sp1other.CopyFrom(sp1);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: sp2other = sp2;
			sp2other.CopyFrom(sp2);
		}

		Vec < 3> v1 = mesh[sp1] - mesh[pi1], v2 = mesh[sp2] - mesh[pi1], v3 = mesh[sp1] - mesh[pi2], v4 = mesh[sp2] - mesh[pi2];
		double vol = 0.5 * (netgen.GlobalMembers.Cross(v1,v2).Length() + netgen.GlobalMembers.Cross(v3,v4).Length());
		h = ngsimd.GlobalMembers.sqrt(vol);
		h = 0;

		sbad = CalcTriangleBadness(mesh[pi1],mesh[pi2],mesh[sp1],0,0) + CalcTriangleBadness(mesh[pi2],mesh[pi1],mesh[sp2],0,0);



		bool puretet = true;
		for (int k = 0; puretet && k < hasbothpoints.Size(); k++)
		{
		  if (mesh[hasbothpoints[k]].GetType() != TET)
		  {
			puretet = false;
		  }
		}
		for (int k = 0; puretet && k < hasbothpointsother.Size(); k++)
		{
		  if (mesh[hasbothpointsother[k]].GetType() != TET)
		  {
			puretet = false;
		  }
		}
		if (!puretet)
		{
		  continue;
		}

		int nsuround = hasbothpoints.Size();
		int nsuroundother = hasbothpointsother.Size();

		Array< int > outerpoints = new Array< int >(nsuround + 1);
		outerpoints[0] = sp1;

		for (int i = 0; i < nsuround; i++)
		{
			bool done = false;
			for (int jj = i; !done && jj < hasbothpoints.Size(); jj++)
			{
			for (int k = 0; !done && k < 4; k++)
			{
			  if (mesh[hasbothpoints[jj]][k] == outerpoints[i])
			  {
			  done = true;
			  for (int l = 0; l < 4; l++)
			  {
				if (mesh[hasbothpoints[jj]][l] != pi1 && mesh[hasbothpoints[jj]][l] != pi2 && mesh[hasbothpoints[jj]][l] != outerpoints[i])
				{
				  outerpoints[i + 1] = mesh[hasbothpoints[jj]][l];
				}
			  }
			  }
			}
			if (done)
			{
				ElementIndex aux = hasbothpoints[i];
				hasbothpoints[i] = hasbothpoints[jj];
				hasbothpoints[jj] = aux;
			}
			}
		}
		if (outerpoints[nsuround] != sp2)
		{
			cerr << "OJE OJE OJE" << "\n";
			(*testout) << "OJE OJE OJE" << "\n";
			(*testout) << "hasbothpoints: " << "\n";
			for (int ii = 0; ii < hasbothpoints.Size(); ii++)
			{
			(*testout) << mesh[hasbothpoints[ii]] << "\n";
			for (int jj = 0; jj < mesh[hasbothpoints[ii]].GetNP(); jj++)
			{
			  if (mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] > 0)
			  {
				(*testout) << mesh[hasbothpoints[ii]][jj] << " between " << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] << " and " << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][1] << "\n";
			  }
			}
			}
			(*testout) << "outerpoints: " << outerpoints << "\n";
			(*testout) << "sel1 " << mesh[sel1] << "\n" << "sel2 " << mesh[sel2] << "\n";
			for (int ii = 0; ii < 3; ii++)
			{
			if (mesh.mlbetweennodes[mesh[sel1][ii]][0] > 0)
			{
			  (*testout) << mesh[sel1][ii] << " between " << mesh.mlbetweennodes[mesh[sel1][ii]][0] << " and " << mesh.mlbetweennodes[mesh[sel1][ii]][1] << "\n";
			}
			if (mesh.mlbetweennodes[mesh[sel2][ii]][0] > 0)
			{
			  (*testout) << mesh[sel2][ii] << " between " << mesh.mlbetweennodes[mesh[sel2][ii]][0] << " and " << mesh.mlbetweennodes[mesh[sel2][ii]][1] << "\n";
			}
			}
		}


		Array< int > outerpointsother = new Array< int >();

		if (nsuroundother > 0)
		{
			outerpointsother.SetSize(nsuroundother + 1);
			outerpointsother[0] = sp2other;
		}

		for (int i = 0; i < nsuroundother; i++)
		{
			bool done = false;
			for (int jj = i; !done && jj < hasbothpointsother.Size(); jj++)
			{
			for (int k = 0; !done && k < 4; k++)
			{
			  if (mesh[hasbothpointsother[jj]][k] == outerpointsother[i])
			  {
			  done = true;
			  for (int l = 0; l < 4; l++)
			  {
				if (mesh[hasbothpointsother[jj]][l] != pi1other && mesh[hasbothpointsother[jj]][l] != pi2other && mesh[hasbothpointsother[jj]][l] != outerpointsother[i])
				{
				  outerpointsother[i + 1] = mesh[hasbothpointsother[jj]][l];
				}
			  }
			  }
			}
			if (done)
			{
				ElementIndex aux = hasbothpointsother[i];
				hasbothpointsother[i] = hasbothpointsother[jj];
				hasbothpointsother[jj] = aux;
			}
			}
		}
		if (nsuroundother > 0 && outerpointsother[nsuroundother] != sp1other)
		{
			cerr << "OJE OJE OJE (other)" << "\n";
			(*testout) << "OJE OJE OJE (other)" << "\n";
			(*testout) << "pi1 " << pi1 << " pi2 " << pi2 << " sp1 " << sp1 << " sp2 " << sp2 << "\n";
			(*testout) << "hasbothpoints: " << "\n";
			for (int ii = 0; ii < hasbothpoints.Size(); ii++)
			{
			(*testout) << mesh[hasbothpoints[ii]] << "\n";
			for (int jj = 0; jj < mesh[hasbothpoints[ii]].GetNP(); jj++)
			{
			  if (mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] > 0)
			  {
				(*testout) << mesh[hasbothpoints[ii]][jj] << " between " << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] << " and " << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][1] << "\n";
			  }
			}
			}
			(*testout) << "outerpoints: " << outerpoints << "\n";
			(*testout) << "sel1 " << mesh[sel1] << "\n" << "sel2 " << mesh[sel2] << "\n";
			for (int ii = 0; ii < 3; ii++)
			{
			if (mesh.mlbetweennodes[mesh[sel1][ii]][0] > 0)
			{
			  (*testout) << mesh[sel1][ii] << " between " << mesh.mlbetweennodes[mesh[sel1][ii]][0] << " and " << mesh.mlbetweennodes[mesh[sel1][ii]][1] << "\n";
			}
			if (mesh.mlbetweennodes[mesh[sel2][ii]][0] > 0)
			{
			  (*testout) << mesh[sel2][ii] << " between " << mesh.mlbetweennodes[mesh[sel2][ii]][0] << " and " << mesh.mlbetweennodes[mesh[sel2][ii]][1] << "\n";
			}
			}

			(*testout) << "pi1other " << pi1other << " pi2other " << pi2other << " sp1other " << sp1other << " sp2other " << sp2other << "\n";
			(*testout) << "hasbothpointsother: " << "\n";
			for (int ii = 0; ii < hasbothpointsother.Size(); ii++)
			{
			(*testout) << mesh[hasbothpointsother[ii]] << "\n";
			for (int jj = 0; jj < mesh[hasbothpointsother[ii]].GetNP(); jj++)
			{
			  if (mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][0] > 0)
			  {
				(*testout) << mesh[hasbothpointsother[ii]][jj] << " between " << mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][0] << " and " << mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][1] << "\n";
			  }
			}
			}
			(*testout) << "outerpoints: " << outerpointsother << "\n";
			(*testout) << "sel1other " << mesh[sel1other] << "\n" << "sel2other " << mesh[sel2other] << "\n";
			for (int ii = 0; ii < 3; ii++)
			{
			if (mesh.mlbetweennodes[mesh[sel1other][ii]][0] > 0)
			{
			  (*testout) << mesh[sel1other][ii] << " between " << mesh.mlbetweennodes[mesh[sel1other][ii]][0] << " and " << mesh.mlbetweennodes[mesh[sel1other][ii]][1] << "\n";
			}
			if (mesh.mlbetweennodes[mesh[sel2other][ii]][0] > 0)
			{
			  (*testout) << mesh[sel2other][ii] << " between " << mesh.mlbetweennodes[mesh[sel2other][ii]][0] << " and " << mesh.mlbetweennodes[mesh[sel2other][ii]][1] << "\n";
			}
			}
		}

		bad1 = 0;
		for (int i = 0; i < hasbothpoints.Size(); i++)
		{
		  bad1 += CalcBad(mesh.Points(), mesh[hasbothpoints[i]],h);
		}
		for (int i = 0; i < hasbothpointsother.Size(); i++)
		{
		  bad1 += CalcBad(mesh.Points(), mesh[hasbothpointsother[i]],h);
		}
		bad1 /= (double)(hasbothpoints.Size() + hasbothpointsother.Size());


		int startpoints;
		int startpointsother;


		if (outerpoints.Size() == 3)
		{
		  startpoints = 1;
		}
		else if (outerpoints.Size() == 4)
		{
		  startpoints = 2;
		}
		else
		{
		  startpoints = outerpoints.Size();
		}

		if (outerpointsother.Size() == 3)
		{
		  startpointsother = 1;
		}
		else if (outerpointsother.Size() == 4)
		{
		  startpointsother = 2;
		}
		else
		{
		  startpointsother = outerpointsother.Size();
		}


		Array< Array < Element >  > newelts = new Array< Array < Element >  >(startpoints);
		Array< Array < Element >  > neweltsother = new Array< Array < Element >  >(startpointsother);

		double minbad = 1e50;
		double minbadother = 1e50;
		double currbad;
		int minpos = -1;
		int minposother = -1;

		//(*testout) << "pi1 " << pi1 << " pi2 " << pi2 << " outerpoints " << outerpoints << endl;

		for (int i = 0; i < startpoints; i++)
		{
			newelts[i] = new Array <Element>(2 * (nsuround - 1));

			for (int jj = 0; jj < nsuround - 1; jj++)
			{
			(*newelts[i])[2 * jj] = new Element(TET);
			(*newelts[i])[2 * jj + 1] = new Element(TET);
			Element newel1 = *((*newelts[i])[2 * jj]);
			Element newel2 = *((*newelts[i])[2 * jj + 1]);

			newel1[0] = pi1;
			newel1[1] = outerpoints[i];
			newel1[2] = outerpoints[(i + jj + 1) % outerpoints.Size()];
			newel1[3] = outerpoints[(i + jj + 2) % outerpoints.Size()];

			newel2[0] = pi2;
			newel2[1] = outerpoints[i];
			newel2[2] = outerpoints[(i + jj + 2) % outerpoints.Size()];
			newel2[3] = outerpoints[(i + jj + 1) % outerpoints.Size()];


			//(*testout) << "j " << j << " newel1 " << newel1[0] << " "<< newel1[1] << " "<< newel1[2] << " "<< newel1[3] << endl
			//     << " newel2 " << newel2[0] << " "<< newel2[1] << " "<< newel2[2] << " "<< newel2[3] << endl;

			newel1.SetIndex(mattype);
			newel2.SetIndex(mattype);

			}

			bool wrongorientation = true;
			for (int jj = 0; wrongorientation && jj < newelts[i].Size(); jj++)
			{
		  wrongorientation = wrongorientation && WrongOrientation(mesh.Points(), *(*newelts[i])[jj]);
			}

			currbad = 0;

			for (int jj = 0; jj < newelts[i].Size(); jj++)
			{
			if (wrongorientation)
			{
			  netgen.GlobalMembers.Swap(ref (*(*newelts[i])[jj])[2], ref (*(*newelts[i])[jj])[3]);
			}


			// not two new faces on same surface
			Array<int> face_index = new Array<int>();
			for (int k = 0; k < surfaceindicesonnode[(*(*newelts[i])[jj])[0]].Size(); k++)
			{
			  face_index.Append(surfaceindicesonnode[(*(*newelts[i])[jj])[0]][k]);
			}

			for (int k = 1; k < 4; k++)
			{
				for (int l = 0; l < face_index.Size(); l++)
				{
				if (face_index[l] != -1 && !(surfaceindicesonnode[(*(*newelts[i])[jj])[k]].Contains(face_index[l])))
				{
				  face_index[l] = -1;
				}
				}

			}

			for (int k = 0; k < face_index.Size(); k++)
			{
			  if (face_index[k] != -1)
			  {
				currbad += 1e12;
			  }
			}


			currbad += CalcBad(mesh.Points(),*(*newelts[i])[jj],h);


			}

			//currbad /= double(newelts[i]->Size());



			if (currbad < minbad)
			{
			minbad = currbad;
			minpos = i;
			}

		}

		if (startpointsother == 0)
		{
		  minbadother = 0;
		}

		for (int i = 0; i < startpointsother; i++)
		{
			neweltsother[i] = new Array <Element>(2 * (nsuroundother));

			for (int jj = 0; jj < nsuroundother; jj++)
			{
			(*neweltsother[i])[2 * jj] = new Element(TET);
			(*neweltsother[i])[2 * jj + 1] = new Element(TET);
			Element newel1 = *((*neweltsother[i])[2 * jj]);
			Element newel2 = *((*neweltsother[i])[2 * jj + 1]);

			newel1[0] = pi1other;
			newel1[1] = outerpointsother[i];
			newel1[2] = outerpointsother[(i + jj + 1) % outerpointsother.Size()];
			newel1[3] = outerpointsother[(i + jj + 2) % outerpointsother.Size()];

			newel2[0] = pi2other;
			newel2[1] = outerpointsother[i];
			newel2[2] = outerpointsother[(i + jj + 2) % outerpointsother.Size()];
			newel2[3] = outerpointsother[(i + jj + 1) % outerpointsother.Size()];


			//(*testout) << "j " << j << " newel1 " << newel1[0] << " "<< newel1[1] << " "<< newel1[2] << " "<< newel1[3] << endl
			//	     << " newel2 " << newel2[0] << " "<< newel2[1] << " "<< newel2[2] << " "<< newel2[3] << endl;

			newel1.SetIndex(othermattype);
			newel2.SetIndex(othermattype);

			}

			bool wrongorientation = true;
			for (int jj = 0; wrongorientation && jj < neweltsother[i].Size(); jj++)
			{
		  wrongorientation = wrongorientation && WrongOrientation(mesh.Points(), *(*neweltsother[i])[jj]);
			}

			currbad = 0;

			for (int jj = 0; jj < neweltsother[i].Size(); jj++)
			{
			if (wrongorientation)
			{
			  netgen.GlobalMembers.Swap(ref (*(*neweltsother[i])[jj])[2], ref (*(*neweltsother[i])[jj])[3]);
			}

			currbad += CalcBad(mesh.Points(),*(*neweltsother[i])[jj],h);
			}

			//currbad /= double(neweltsother[i]->Size());



			if (currbad < minbadother)
			{
			minbadother = currbad;
			minposother = i;
			}

		}

		//(*testout) << "minbad " << minbad << " bad1 " << bad1 << endl;


		double sbadnew = CalcTriangleBadness(mesh[pi1],mesh[sp2],mesh[sp1],0,0) + CalcTriangleBadness(mesh[pi2],mesh[sp1],mesh[sp2],0,0);


		int denom = newelts[minpos].Size();
		if (minposother >= 0)
		{
		  denom += neweltsother[minposother].Size();
		}


		if ((minbad + minbadother) / (double)denom < bad1 && sbadnew < sbad)
		{
			cnt++;

			swapped = true;


			int start1 = -1;
			for (int l = 0; l < 3; l++)
			{
		  if (mesh[sel1][l] == pi1)
		  {
			start1 = l;
		  }
			}
			if (mesh[sel1][(start1 + 1) % 3] == pi2)
			{
			mesh[sel1][0] = pi1;
			mesh[sel1][1] = sp2;
			mesh[sel1][2] = sp1;
			mesh[sel2][0] = pi2;
			mesh[sel2][1] = sp1;
			mesh[sel2][2] = sp2;
			}
			else
			{
			mesh[sel1][0] = pi2;
			mesh[sel1][1] = sp2;
			mesh[sel1][2] = sp1;
			mesh[sel2][0] = pi1;
			mesh[sel2][1] = sp1;
			mesh[sel2][2] = sp2;
			}
			//(*testout) << "changed surface element " << sel1 << " to " << mesh[sel1] << ", " << sel2 << " to " << mesh[sel2] << endl;

			for (int l = 0; l < 3; l++)
			{
			surfaceelementsonnode.Add(mesh[sel1][l],sel1);
			surfaceelementsonnode.Add(mesh[sel2][l],sel2);
			}



			if (periodic)
			{
			start1 = -1;
			for (int l = 0; l < 3; l++)
			{
			  if (mesh[sel1other][l] == pi1other)
			  {
				start1 = l;
			  }
			}



			//(*testout) << "changed surface elements " << mesh[sel1other] << " and " << mesh[sel2other] << endl;
			if (mesh[sel1other][(start1 + 1) % 3] == pi2other)
			{
				mesh[sel1other][0] = pi1other;
				mesh[sel1other][1] = sp2other;
				mesh[sel1other][2] = sp1other;
				mesh[sel2other][0] = pi2other;
				mesh[sel2other][1] = sp1other;
				mesh[sel2other][2] = sp2other;
				//(*testout) << "       with rule 1" << endl;
			}
			else
			{
				mesh[sel1other][0] = pi2other;
				mesh[sel1other][1] = sp2other;
				mesh[sel1other][2] = sp1other;
				mesh[sel2other][0] = pi1other;
				mesh[sel2other][1] = sp1other;
				mesh[sel2other][2] = sp2other;
				//(*testout) << "       with rule 2" << endl;
			}
			//(*testout) << "         to " << mesh[sel1other] << " and " << mesh[sel2other] << endl;

			//(*testout) << "  and surface element " << sel1other << " to " << mesh[sel1other] << ", " << sel2other << " to " << mesh[sel2other] << endl;

			for (int l = 0; l < 3; l++)
			{
				surfaceelementsonnode.Add(mesh[sel1other][l],sel1other);
				surfaceelementsonnode.Add(mesh[sel2other][l],sel2other);
			}
			}




			for (int i = 0; i < hasbothpoints.Size(); i++)
			{
			mesh[hasbothpoints[i]] = *(*newelts[minpos])[i];

			for (int l = 0; l < 4; l++)
			{
			  elementsonnode.Add((*(*newelts[minpos])[i])[l],hasbothpoints[i]);
			}
			}

			for (int i = hasbothpoints.Size(); i < (*newelts[minpos]).Size(); i++)
			{
			ElementIndex ni = mesh.AddVolumeElement(*(*newelts[minpos])[i]);

			for (int l = 0; l < 4; l++)
			{
			  elementsonnode.Add((*(*newelts[minpos])[i])[l],ni);
			}
			}

			if (hasbothpointsother.Size() > 0)
			{
			for (int i = 0; i < hasbothpointsother.Size(); i++)
			{
				mesh[hasbothpointsother[i]] = *(*neweltsother[minposother])[i];
				for (int l = 0; l < 4; l++)
				{
			  elementsonnode.Add((*(*neweltsother[minposother])[i])[l],hasbothpointsother[i]);
				}
			}

			for (int i = hasbothpointsother.Size(); i < (*neweltsother[minposother]).Size(); i++)
			{
				ElementIndex ni = mesh.AddVolumeElement(*(*neweltsother[minposother])[i]);
				for (int l = 0; l < 4; l++)
				{
			  elementsonnode.Add((*(*neweltsother[minposother])[i])[l],ni);
				}
			}
			}



		}

		for (int i = 0; i < newelts.Size(); i++)
		{
			for (int jj = 0; jj < newelts[i].Size(); jj++)
			{
		  newelts[i][jj] = null;
			}
			newelts[i] = null;
		}

		for (int i = 0; i < neweltsother.Size(); i++)
		{
			for (int jj = 0; jj < neweltsother[i].Size(); jj++)
			{
		  neweltsother[i][jj] = null;
			}
			neweltsother[i] = null;
		}

		}
	}

	PrintMessage(5, cnt, " swaps performed");


	for (int i = 0; i < locidmaps.Size(); i++)
	{
	  locidmaps[i] = null;
	}


	mesh.Compress();

	multithread.task = savetask;
  }


  /*
    2 -> 3 conversion
  */


//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer SwapImprove2_t("MeshOptimize3d::SwapImprove2");

  public void SwapImprove2(Mesh mesh, OPTIMIZEGOAL goal = OPT_QUALITY)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("MeshOptimize3d::SwapImprove2");
	RegionTimer reg = new RegionTimer(SwapImprove2_t);

	PointIndex pi1 = new PointIndex(0);
	PointIndex pi2 = new PointIndex(0);
	PointIndex pi3 = new PointIndex(0);
	PointIndex pi4 = new PointIndex(0);
	PointIndex pi5 = new PointIndex(0);
	Element el21 = new Element(TET);
	Element el22 = new Element(TET);
	Element el31 = new Element(TET);
	Element el32 = new Element(TET);
	Element el33 = new Element(TET);

	int cnt = 0;
	double bad1;
	double bad2;

	int np = mesh.GetNP();
	int ne = mesh.GetNE();
	int nse = mesh.GetNSE();

	if (goal == OPT_CONFORM)
	{
		return;
	}

	// contains at least all elements at node
	TABLE<ElementIndex, PointIndex.BASE> elementsonnode = new TABLE<ElementIndex, PointIndex.BASE>(np);
	TABLE<SurfaceElementIndex, PointIndex.BASE> belementsonnode = new TABLE<SurfaceElementIndex, PointIndex.BASE>(np);

	PrintMessage(3, "SwapImprove2 ");
	(*testout) << "\n" << "Start SwapImprove2" << "\n";
	//  TestOk();


	/*
	  CalcSurfacesOfNode ();
	  for (i = 1; i <= GetNE(); i++)
	  if (volelements.Get(i)[0])
	  if (!mesh.LegalTet (volelements.Get(i)))
	  {
	  cout << "detected illegal tet, 1" << endl;
	  (*testout) << "detected illegal tet1: " << i << endl;
	  }
	*/


	// Calculate total badness

	bad1 = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
	(*testout) << "Total badness = " << bad1 << "\n";
	//  cout << "tot bad = " << bad1 << endl;

	// find elements on node

	for (ElementIndex ei = 0; ei < ne; ei++)
	{
	  for (int j = 0; j < mesh[ei].GetNP(); j++)
	  {
		elementsonnode.Add(mesh[ei][j], ei);
	  }
	}

	for (SurfaceElementIndex sei = 0; sei < nse; sei++)
	{
	  for (int j = 0; j < 3; j++)
	  {
		belementsonnode.Add(mesh[sei][j], sei);
	  }
	}

	for (ElementIndex eli1 = 0; eli1 < ne; eli1++)
	{
		if (multithread.terminate)
		{
	  break;
		}

		if (mesh.ElementType(eli1) == FIXEDELEMENT)
		{
	  continue;
		}

		if (mesh[eli1].GetType() != TET)
		{
	  continue;
		}

		if ((goal == OPT_LEGAL) && mesh.LegalTet(mesh[eli1]) && CalcBad(mesh.Points(), mesh[eli1], 0) < 1e3)
		{
	  continue;
		}

		if (mesh.GetDimension() == 3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(eli1).GetIndex())
		{
			continue;
		}

		// cout << "eli = " << eli1 << endl;
		//      (*testout) << "swapimp2, eli = " << eli1 << "; el = " << mesh[eli1] << endl;

		for (int j = 0; j < 4; j++)
		{
		// loop over faces

		Element elem = mesh[eli1];
		// if (elem[0] < PointIndex::BASE) continue;
		if (elem.IsDeleted())
		{
			continue;
		}

		int mattyp = elem.GetIndex();

		switch (j)
		{
		  case 0:
			pi1 = elem.PNum(1);
			pi2 = elem.PNum(2);
			pi3 = elem.PNum(3);
			pi4 = elem.PNum(4);
			break;
		  case 1:
			pi1 = elem.PNum(1);
			pi2 = elem.PNum(4);
			pi3 = elem.PNum(2);
			pi4 = elem.PNum(3);
			break;
		  case 2:
			pi1 = elem.PNum(1);
			pi2 = elem.PNum(3);
			pi3 = elem.PNum(4);
			pi4 = elem.PNum(2);
			break;
		  case 3:
			pi1 = elem.PNum(2);
			pi2 = elem.PNum(4);
			pi3 = elem.PNum(3);
			pi4 = elem.PNum(1);
			break;
		}


		bool bface = false;
		for (int k = 0; k < belementsonnode[pi1].Size(); k++)
		{
			Element2d bel = mesh[belementsonnode[pi1][k]];

			bool bface1 = true;
			for (int l = 0; l < 3; l++)
			{
		  if (bel[l] != pi1 && bel[l] != pi2 && bel[l] != pi3)
		  {
			  bface1 = false;
			  break;
		  }
			}

			if (bface1)
			{
			bface = true;
			break;
			}
		}

		if (bface)
		{
			continue;
		}


		FlatArray<ElementIndex> row = elementsonnode[pi1];
		for (int k = 0; k < row.Size(); k++)
		{
			ElementIndex eli2 = row[k];

			// cout << "\rei1 = " << eli1 << ", pi1 = " << pi1 << ", k = " << k << ", ei2 = " << eli2
			// << ", getne = " << mesh.GetNE();

			if (eli1 != eli2)
			{
			Element elem2 = mesh[eli2];
			if (elem2.IsDeleted())
			{
				continue;
			}
			if (elem2.GetType() != TET)
			{
			  continue;
			}

			int comnodes = 0;
			for (int l = 1; l <= 4; l++)
			{
			  if (elem2.PNum(l) == pi1 || elem2.PNum(l) == pi2 || elem2.PNum(l) == pi3)
			  {
			  comnodes++;
			  }
			  else
			  {
			  pi5 = elem2.PNum(l);
			  }
			}

			if (comnodes == 3)
			{
				bad1 = CalcBad(mesh.Points(), elem, 0) + CalcBad(mesh.Points(), elem2, 0);

				if (!mesh.LegalTet(elem) || !mesh.LegalTet(elem2))
				{
			  bad1 += 1e4;
				}


				el31.PNum(1) = pi1;
				el31.PNum(2) = pi2;
				el31.PNum(3) = pi5;
				el31.PNum(4) = pi4;
				el31.SetIndex(mattyp);

				el32.PNum(1) = pi2;
				el32.PNum(2) = pi3;
				el32.PNum(3) = pi5;
				el32.PNum(4) = pi4;
				el32.SetIndex(mattyp);

				el33.PNum(1) = pi3;
				el33.PNum(2) = pi1;
				el33.PNum(3) = pi5;
				el33.PNum(4) = pi4;
				el33.SetIndex(mattyp);

				bad2 = CalcBad(mesh.Points(), el31, 0) + CalcBad(mesh.Points(), el32, 0) + CalcBad(mesh.Points(), el33, 0);


				el31.flags.illegal_valid = 0;
				el32.flags.illegal_valid = 0;
				el33.flags.illegal_valid = 0;

				if (!mesh.LegalTet(el31) || !mesh.LegalTet(el32) || !mesh.LegalTet(el33))
				{
			  bad2 += 1e4;
				}


				bool do_swap = (bad2 < bad1);

				if (((bad2 < 1e6) || (bad2 < 10 * bad1)) && mesh.BoundaryEdge(pi4, pi5))
				{
			  do_swap = true;
				}

				if (do_swap)
				{
				//			  cout << "do swap, eli1 = " << eli1 << "; eli2 = " << eli2 << endl;
				//			  (*mycout) << "2->3 " << flush;
				cnt++;

				el31.flags.illegal_valid = 0;
				el32.flags.illegal_valid = 0;
				el33.flags.illegal_valid = 0;

				mesh[eli1] = el31;
				mesh[eli2] = el32;

				ElementIndex neli = mesh.AddVolumeElement(el33);

				/*
				  if (!LegalTet(el31) || !LegalTet(el32) ||
				  !LegalTet(el33))
				  {
				  cout << "Swap to illegal tets !!!" << endl;
				  }
				*/
				// cout << "neli = " << neli << endl;
				for (int l = 0; l < 4; l++)
				{
					elementsonnode.Add(el31[l], eli1);
					elementsonnode.Add(el32[l], eli2);
					elementsonnode.Add(el33[l], neli);
				}

				break;
				}
			}
			}
		}
		}
	}


	PrintMessage(5, cnt, " swaps performed");



	/*
	  CalcSurfacesOfNode ();
	  for (i = 1; i <= GetNE(); i++)
	  if (volelements.Get(i).PNum(1))
	  if (!LegalTet (volelements.Get(i)))
	  {
	  cout << "detected illegal tet, 2" << endl;
	  (*testout) << "detected illegal tet2: " << i << endl;
	  }
	*/


	bad1 = CalcTotalBad(mesh.Points(), mesh.VolumeElements());
	(*testout) << "Total badness = " << bad1 << "\n";
	(*testout) << "swapimprove2 done" << "\n";
	//  (*mycout) << "Vol = " << CalcVolume (points, volelements) << "\n";
  }

  public double CalcBad(Mesh.T_POINTS points, Element elem, double h)
  {
	if (elem.GetType() == TET)
	{
	  return CalcTetBadness(points[elem[0]], points[elem[1]], points[elem[2]], points[elem[3]], h, mp);
	}
	return 0;
  }


  public double CalcTotalBad(Mesh.T_POINTS points, Array<Element, 0, uint> elements)
  {
	return netgen.CalcTotalBad(points, elements, mp);
  }

}


/* Functional depending of inner point inside triangular surface */


public class MinFunctionSum : MinFunction
{
  protected Array<MinFunction> functions = new Array<MinFunction>();


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double Func(Vector x);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void Grad(const Vector & x, Vector & g) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual void Grad(Vector x, Vector g);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & g) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double FuncGrad(Vector x, Vector g);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double FuncDeriv(Vector x, Vector dir, ref double deriv);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double GradStopping(const Vector & x) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double GradStopping(Vector x);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void AddFunction(MinFunction fun);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MinFunction & Function(int i) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  MinFunction Function(int i);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  MinFunction Function(int i);
}


public class PointFunction1 : MinFunction
{
  private Mesh.T_POINTS points;
  private readonly Array<INDEX_3> faces;
  private readonly MeshingParameters mp;
  private double h;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  PointFunction1(Mesh::T_POINTS apoints, Array<INDEX_3> afaces, MeshingParameters amp, double ah);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double Func(Vector x);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double FuncDeriv(Vector x, Vector dir, ref double deriv);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & g) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double FuncGrad(Vector x, Vector g);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double GradStopping(const Vector & x) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double GradStopping(Vector x);
}

public class JacobianPointFunction : MinFunction
{
  public Mesh.T_POINTS points;
  public readonly Array<Element, 0, uint> elements;
  public TABLE<INDEX> elementsonpoint = new TABLE<INDEX>();
  public PointIndex actpind = new PointIndex();

  public bool onplane;
  public Vec < 3> nv;

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  JacobianPointFunction(Mesh::T_POINTS apoints, Array<Element, 0, uint> aelements);
  public virtual void Dispose()
  {
	  ;
  }
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual void SetPointIndex(PointIndex aactpind);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double Func(Vector x);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & g) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double FuncGrad(Vector x, Vector g);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual double FuncDeriv(Vector x, Vector dir, ref double deriv);

  public void SetNV(Vec < 3> anv)
  {
	  nv = anv;
	  onplane = true;
  }
  public void UnSetNV()
  {
	  onplane = false;
  }
}





#if SOLIDGEOM
#endif

namespace netgen
{


/*
  void Mesh :: SwapImprove2 (OPTIMIZEGOAL goal)
  {
  int i, j;
  int eli1, eli2;
  int mattyp;

  Element el31(4), el32(4), el33(4);
  double bad1, bad2;


  INDEX_3_HASHTABLE<INDEX_2> elsonface (GetNE());

  (*mycout) << "SwapImprove2 " << endl;
  (*testout) << "\n" << "Start SwapImprove2" << "\n";

  // Calculate total badness

  if (goal == OPT_QUALITY)
  {
  double bad1 = CalcTotalBad (points, volelements);
  (*testout) << "Total badness = " << bad1 << endl;
  }

  // find elements on node


  Element2d face;
  for (i = 1; i <= GetNE(); i++)
  if ( (i > eltyps.Size()) || (eltyps.Get(i) != FIXEDELEMENT) )
  {
  const Element & el = VolumeElement(i);
  if (!el.PNum(1)) continue;

  for (j = 1; j <= 4; j++)
  {
  el.GetFace (j, face);
  INDEX_3 i3 (face.PNum(1), face.PNum(2), face.PNum(3));
  i3.Sort();


  int bnr, posnr;
  if (!elsonface.PositionCreate (i3, bnr, posnr))
  {
  INDEX_2 i2;
  elsonface.GetData (bnr, posnr, i3, i2);
  i2.I2() = i;
  elsonface.SetData (bnr, posnr, i3, i2);
  }
  else
  {
  INDEX_2 i2 (i, 0);
  elsonface.SetData (bnr, posnr, i3, i2);
  }

  //  	    if (elsonface.Used (i3))
  //  	      {
  //  		INDEX_2 i2 = elsonface.Get(i3);
  //  		i2.I2() = i;
  //  		elsonface.Set (i3, i2);
  //  	      }
  //  	    else
  //  	      {
  //  		INDEX_2 i2 (i, 0);
  //  		elsonface.Set (i3, i2);
  //  	      }

  }
  }

  BitArray original(GetNE());
  original.Set();

  for (i = 1; i <= GetNSE(); i++)
  {
  const Element2d & sface = SurfaceElement(i);
  INDEX_3 i3 (sface.PNum(1), sface.PNum(2), sface.PNum(3));
  i3.Sort();
  INDEX_2 i2(0,0);
  elsonface.Set (i3, i2);
  }


  for (i = 1; i <= elsonface.GetNBags(); i++)
  for (j = 1; j <= elsonface.GetBagSize(i); j++)
  {
  INDEX_3 i3;
  INDEX_2 i2;
  elsonface.GetData (i, j, i3, i2);


  int eli1 = i2.I1();
  int eli2 = i2.I2();

  if (eli1 && eli2 && original.Test(eli1) && original.Test(eli2) )
  {
  Element & elem = volelements.Elem(eli1);
  Element & elem2 = volelements.Elem(eli2);

  int pi1 = i3.I1();
  int pi2 = i3.I2();
  int pi3 = i3.I3();

  int pi4 = elem.PNum(1) + elem.PNum(2) + elem.PNum(3) + elem.PNum(4) - pi1 - pi2 - pi3;
  int pi5 = elem2.PNum(1) + elem2.PNum(2) + elem2.PNum(3) + elem2.PNum(4) - pi1 - pi2 - pi3;






  el31.PNum(1) = pi1;
  el31.PNum(2) = pi2;
  el31.PNum(3) = pi3;
  el31.PNum(4) = pi4;
  el31.SetIndex (mattyp);
	    
  if (WrongOrientation (points, el31))
  swap (pi1, pi2);


  bad1 = CalcBad (points, elem, 0) + 
  CalcBad (points, elem2, 0); 
	    
  //	    if (!LegalTet(elem) || !LegalTet(elem2))
  //	      bad1 += 1e4;

	    
  el31.PNum(1) = pi1;
  el31.PNum(2) = pi2;
  el31.PNum(3) = pi5;
  el31.PNum(4) = pi4;
  el31.SetIndex (mattyp);
	    
  el32.PNum(1) = pi2;
  el32.PNum(2) = pi3;
  el32.PNum(3) = pi5;
  el32.PNum(4) = pi4;
  el32.SetIndex (mattyp);
		      
  el33.PNum(1) = pi3;
  el33.PNum(2) = pi1;
  el33.PNum(3) = pi5;
  el33.PNum(4) = pi4;
  el33.SetIndex (mattyp);
	    
  bad2 = CalcBad (points, el31, 0) + 
  CalcBad (points, el32, 0) +
  CalcBad (points, el33, 0); 
	    
  //	    if (!LegalTet(el31) || !LegalTet(el32) ||
  //		!LegalTet(el33))
  //	      bad2 += 1e4;
	    
	    
  int swap = (bad2 < bad1);

  INDEX_2 hi2b(pi4, pi5);
  hi2b.Sort();
	    
  if ( ((bad2 < 1e6) || (bad2 < 10 * bad1)) &&
  boundaryedges->Used (hi2b) )
  swap = 1;
	    
  if (swap)
  {
  (*mycout) << "2->3 " << flush;
		
  volelements.Elem(eli1) = el31;
  volelements.Elem(eli2) = el32;
  volelements.Append (el33);
		
  original.Clear (eli1);
  original.Clear (eli2);
  }
  }
  }
  
  (*mycout) << endl;

  if (goal == OPT_QUALITY)
  {
  bad1 = CalcTotalBad (points, volelements);
  (*testout) << "Total badness = " << bad1 << endl;
  }

  //  FindOpenElements ();

  (*testout) << "swapimprove2 done" << "\n";
  }

*/
}
