namespace netgen
{

	public class Meshing2
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void BlockFillLocalH(Mesh mesh, MeshingParameters mp)
		  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("Meshing2::BlockFill");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1 = NgProfiler::CreateTimer("Meshing2::BlockFill 1");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer2 = NgProfiler::CreateTimer("Meshing2::BlockFill 2");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer3 = NgProfiler::CreateTimer("Meshing2::BlockFill 3");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer4 = NgProfiler::CreateTimer("Meshing2::BlockFill 4");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(BlockFillLocalH_timer);
        
			NgProfiler.StartTimer(BlockFillLocalH_timer1);
        
			double filldist = mp.filldist;
        
			Console.Write("blockfill local h");
			Console.Write("\n");
			Console.Write("rel filldist = ");
			Console.Write(filldist);
			Console.Write("\n");
			PrintMessage(3, "blockfill local h");
        
			Array<Point < 3>> npoints = new Array<Point < 3>>();
        
			// adfront -> CreateTrees();
        
			Box < 3> bbox(Box < 3>.EMPTY_BOX);
			double maxh = 0;
        
			for (int i = 0; i < adfront.GetNFL(); i++)
			{
			FrontLine line = adfront.GetLine(i);
        
			const Point < 3> & p1 = adfront.GetPoint(line.L().I1());
			const Point < 3> & p2 = adfront.GetPoint(line.L().I2());
        
				maxh = Math.Max(maxh, Dist(p1, p2));
        
			bbox.Add(p1);
			bbox.Add(p2);
			}
        
        
			Console.Write("bbox = ");
			Console.Write(bbox);
			Console.Write("\n");
        
        
			// Point<3> mpc = bbox.Center();
			bbox.Increase(bbox.Diam() / 2);
			Box < 3> meshbox = bbox;
        
			NgProfiler.StopTimer(BlockFillLocalH_timer1);
			NgProfiler.StartTimer(BlockFillLocalH_timer2);
        
        
			LocalH loch2 = new LocalH(bbox, 1, 2);
        
			if (mp.maxh < maxh)
			{
				maxh = mp.maxh;
			}
        
			bool changed;
			do
			{
			mesh.LocalHFunction().ClearFlags();
        
			for (int i = 0; i < adfront.GetNFL(); i++)
			{
				FrontLine line = adfront.GetLine(i);
        
				Box < 3> bbox(adfront.GetPoint(line.L().I1()));
				bbox.Add(adfront.GetPoint(line.L().I2()));
        
        
				double filld = filldist * bbox.Diam();
				bbox.Increase(filld);
        
				mesh.LocalHFunction().CutBoundary(bbox);
			}
        
        
			mesh.LocalHFunction().FindInnerBoxes(adfront, null);
        
			npoints.SetSize(0);
			mesh.LocalHFunction().GetInnerPoints(npoints);
        
			changed = false;
			for (int i = 0; i < npoints.Size(); i++)
			{
				if (mesh.LocalHFunction().GetH(npoints[i]) > 1.2 * maxh)
				{
				mesh.LocalHFunction().SetH(npoints[i], maxh);
				changed = true;
				}
			}
			} while (changed);
        
			NgProfiler.StopTimer(BlockFillLocalH_timer2);
			NgProfiler.StartTimer(BlockFillLocalH_timer3);
        
        
			if (debugparam.slowchecks)
			{
			  (*testout) << "Blockfill with points: " << "\n";
			}
			*testout << "loch = " << mesh.LocalHFunction() << "\n";
        
			*testout << "npoints = " << "\n" << npoints << "\n";
        
			for (int i = 1; i <= npoints.Size(); i++)
			{
			if (meshbox.IsIn(npoints.Get(i)))
			{
				PointIndex gpnum = mesh.AddPoint(npoints.Get(i));
				adfront.AddPoint(npoints.Get(i), gpnum);
        
				if (debugparam.slowchecks)
				{
				(*testout) << npoints.Get(i) << "\n";
        
				Point < 2> p2d(npoints.Get(i)(0), npoints.Get(i)(1));
				if (!adfront.Inside(p2d))
				{
					Console.Write("add outside point");
					Console.Write("\n");
					(*testout) << "outside" << "\n";
				}
				}
        
			}
			}
        
			NgProfiler.StopTimer(BlockFillLocalH_timer3);
			NgProfiler.StartTimer(BlockFillLocalH_timer4);
        
        
		  // find outer points
        
			loch2.ClearFlags();
        
			for (int i = 0; i < adfront.GetNFL(); i++)
			{
			FrontLine line = adfront.GetLine(i);
        
			Box < 3> bbox(adfront.GetPoint(line.L().I1()));
			bbox.Add(adfront.GetPoint(line.L().I2()));
        
			loch2.SetH(bbox.Center(), bbox.Diam());
			}
        
        
			for (int i = 0; i < adfront.GetNFL(); i++)
			{
			FrontLine line = adfront.GetLine(i);
        
			Box < 3> bbox(adfront.GetPoint(line.L().I1()));
			bbox.Add(adfront.GetPoint(line.L().I2()));
        
			bbox.Increase(filldist * bbox.Diam());
			loch2.CutBoundary(bbox);
			}
        
			loch2.FindInnerBoxes(adfront, null);
        
			  // outer points : smooth mesh-grading
			npoints.SetSize(0);
			loch2.GetOuterPoints(npoints);
        
			for (int i = 1; i <= npoints.Size(); i++)
			{
			if (meshbox.IsIn(npoints.Get(i)))
			{
				PointIndex gpnum = mesh.AddPoint(npoints.Get(i));
				adfront.AddPoint(npoints.Get(i), gpnum);
			}
			}
        
			NgProfiler.StopTimer(BlockFillLocalH_timer4);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Delaunay(Mesh mesh, int domainnr, MeshingParameters mp)
		  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("Meshing2::Delaunay - total");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timerstart = NgProfiler::CreateTimer("Meshing2::Delaunay - start");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timerfinish = NgProfiler::CreateTimer("Meshing2::Delaunay - finish");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1 = NgProfiler::CreateTimer("Meshing2::Delaunay - incremental");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1a = NgProfiler::CreateTimer("Meshing2::Delaunay - incremental a");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1b = NgProfiler::CreateTimer("Meshing2::Delaunay - incremental b");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1c = NgProfiler::CreateTimer("Meshing2::Delaunay - incremental c");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1d = NgProfiler::CreateTimer("Meshing2::Delaunay - incremental d");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(Delaunay_timer);
        
        
        
			Console.Write("2D Delaunay meshing (in progress)");
			Console.Write("\n");
        
        
			BlockFillLocalH(mesh, mp);
        
			NgProfiler.StartTimer(Delaunay_timerstart);
        
			// do the delaunay
        
        
			// face bounding box:
			Box < 3> bbox(Box < 3>.EMPTY_BOX);
        
			for (int i = 0; i < adfront.GetNFL(); i++)
			{
			FrontLine line = adfront.GetLine(i);
				bbox.Add(Point < 3> (adfront.GetPoint(line.L [0])));
				bbox.Add(Point < 3> (adfront.GetPoint(line.L [1])));
			}
        
			for (int i = 0; i < mesh.LockedPoints().Size(); i++)
			{
			  bbox.Add(new mesh.Point(mesh.LockedPoints [i]));
			}
        
			Console.Write("bbox = ");
			Console.Write(bbox);
			Console.Write("\n");
        
			// external point
			Vec < 3> vdiag = bbox.PMax() - bbox.PMin();
        
			var old_points = mesh.Points().Range();
			DelaunayTrig startel = new DelaunayTrig();
			startel[0] = mesh.AddPoint(bbox.PMin() + Vec < 3> (-8 * vdiag(0), -8 * vdiag(1), 0));
			startel[1] = mesh.AddPoint(bbox.PMin() + Vec < 3> (+8 * vdiag(0), -8 * vdiag(1), 0));
			startel[2] = mesh.AddPoint(bbox.PMin() + Vec < 3> (0, 8 * vdiag(1), 0));
        
			Box < 3> hbox;
			hbox.Set(mesh[startel[0]]);
			hbox.Add(mesh[startel[1]]);
			hbox.Add(mesh[startel[2]]);
			Point < 3> hp = mesh[startel[0]];
			hp(2) = 1;
			hbox.Add(hp);
			hp(2) = -1;
			hbox.Add(hp);
			BoxTree < 3> searchtree(hbox);
        
			Array<DelaunayTrig> tempels = new Array<DelaunayTrig>();
			startel.CalcCenter(mesh);
        
			tempels.Append(startel);
			searchtree.Insert(startel.BoundingBox(), 0);
        
			Array<int> closeels = new Array<int>();
			Array<int> intersecting = new Array<int>();
			Array<INDEX_2> edges = new Array<INDEX_2>();
        
        
        
        
			// reorder points
			Array<PointIndex, PointIndex.BASE, PointIndex> mixed = new Array<PointIndex, PointIndex.BASE, PointIndex>(old_points.Size());
			int[] prims = {11, 13, 17, 19, 23, 29, 31, 37};
			int prim;
        
			{
			  int i = 0;
			  while (old_points.Size() % prims[i] == 0)
			  {
				  i++;
			  }
			  prim = prims[i];
			}
        
			foreach (PointIndex pi in old_points)
			{
			  mixed[pi] = new PointIndex((prim * pi) % old_points.Size() + PointIndex.BASE);
			}
        
			NgProfiler.StopTimer(Delaunay_timerstart);
			NgProfiler.StartTimer(Delaunay_timer1);
        
        
			foreach (PointIndex i1 in old_points)
			{
				PointIndex i = mixed[i1];
        
				NgProfiler.StartTimer(Delaunay_timer1a);
				Point < 3> newp = mesh[i];
				intersecting.SetSize(0);
				edges.SetSize(0);
        
				searchtree.GetIntersecting(newp, newp, closeels);
				// for (int jj = 0; jj < closeels.Size(); jj++)
				// for (int j = 0; j < tempels.Size(); j++)
				foreach (int j in closeels)
				{
					if (tempels[j][0] < 0)
					{
						continue;
					}
					Point < 3> c = tempels[j].Center();
					double r2 = tempels[j].Radius2();
        
					bool inside = Dist2(mesh[i], c) < r2;
					if (inside)
					{
						intersecting.Append(j);
					}
				}
        
				NgProfiler.StopTimer(Delaunay_timer1a);
				NgProfiler.StartTimer(Delaunay_timer1b);
        
				// find outer edges
				foreach (var j in intersecting)
				{
					DelaunayTrig trig = tempels[j];
					for (int k = 0; k < 3; k++)
					{
						int p1 = trig[k];
						int p2 = trig[(k + 1) % 3];
						INDEX_2 edge = new INDEX_2(p1, p2);
						edge.Sort();
						bool found = false;
						for (int l = 0; l < edges.Size(); l++)
						{
						  if (edges[l] == edge)
						  {
							  edges.Delete(l);
							  found = true;
							  break;
						  }
						}
						if (!found)
						{
							edges.Append(edge);
						}
					}
				}
        
				NgProfiler.StopTimer(Delaunay_timer1b);
				NgProfiler.StartTimer(Delaunay_timer1c);
        
				/*
				for (int j = intersecting.Size()-1; j >= 0; j--)
				  tempels.Delete (intersecting[j]);
				*/
				foreach (int j in intersecting)
				{
					searchtree.DeleteElement(j);
					tempels[j][0] = -1;
					tempels[j][1] = -1;
					tempels[j][2] = -1;
				}
        
				NgProfiler.StopTimer(Delaunay_timer1c);
				NgProfiler.StartTimer(Delaunay_timer1d);
        
				foreach (var edge in edges)
				{
					DelaunayTrig trig = new DelaunayTrig(edge[0], edge[1], i);
					trig.CalcCenter(mesh);
					tempels.Append(trig);
					searchtree.Insert(trig.BoundingBox(), tempels.Size() - 1);
				}
        
				NgProfiler.StopTimer(Delaunay_timer1d);
			}
        
			NgProfiler.StopTimer(Delaunay_timer1);
			NgProfiler.StartTimer(Delaunay_timerfinish);
        
			foreach (DelaunayTrig trig in tempels)
			{
				if (trig[0] < 0)
				{
					continue;
				}
        
				Point < 3> c = Center(mesh[trig[0]], mesh[trig[1]], mesh[trig[2]]);
				if (!adfront.Inside(Point < 2> (c(0),c(1))))
				{
					continue;
				}
        
				Vec < 3> n = Cross(mesh[trig[1]] - mesh[trig[0]], mesh[trig[2]] - mesh[trig[0]]);
				if (n(2) < 0)
				{
					Swap(ref trig[1], ref trig[2]);
				}
        
				Element2d el = new Element2d(trig[0], trig[1], trig[2]);
				el.SetIndex(domainnr);
				mesh.AddSurfaceElement(el);
			}
        
			foreach (PointIndex pi in mesh.Points().Range())
			{
			  *testout << pi << ": " << mesh[pi].Type() << "\n";
			}
        
			NgProfiler.StopTimer(Delaunay_timerfinish);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void AddBoundaryElement(int i1, int i2, PointGeomInfo gi1, PointGeomInfo gi2)
		  {
			//    (*testout) << "add line " << i1 << " - " << i2 << endl;
			if (gi1.trignum == 0 || gi2.trignum == 0)
			{
			PrintSysError("addboundaryelement: illegal geominfo");
			}
			adfront.AddLine(i1 - 1, i2 - 1, gi1, gi2);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void LoadRules(string filename, bool quad)
		{
		  string buf = new string(new char[256]);
		  istream ist;
		  //char *tr1 = NULL;
		  string tr1;
        
		  /*
		  ifstream ist (filename);
		  if (!ist.good())
		    {
		      cerr << "Rule description file " << filename << " not found" << endl;
		      exit (1);
		    }
		  */
        
        
		  if (filename)
		  {
			  //      (*mycout) << "rule-filename = " << filename << endl;
			  ist = new ifstream(filename);
		  }
		  else
		  {
			  /* connect tetrules to one string */
			  string[] hcp;
        
			  // if (!mparam.quad)
			  if (!quad)
			  {
			  hcp = triarules;
			  PrintMessage(3, "load internal triangle rules");
			  }
			  else
			  {
			  hcp = quadrules;
			  PrintMessage(3, "load internal quad rules");
			  // LoadRules ("rules/quad.rls");
			  }
        
			  uint len = 0;
			  while hcp
			  {
			  //	  (*testout) << "POS2 *hcp " << *hcp << endl;
			  len += Convert.ToStringhcp.Length;
			  hcp++;
			  }
			  //tr1 = new char[len+1];
			  //tr1[0] = 0;
			  tr1.reserve(len + 1);
        
        
			  // if (!mparam.quad)
			  if (!quad)
			  {
			hcp = triarules;
			  }
			  else
			  {
			hcp = quadrules;
			  }
        
        
			  //char * tt1 = tr1;
			  while hcp
			  {
			  //strcat (tt1, *hcp);
			  //tt1 += strlen (*hcp);
			  tr1.appendhcp;
			  hcp++;
			  }
        
		#if WIN32
			  // VC++ 2005 workaround
			  for (int i = 0; i < tr1.Length; i++)
			  {
			if (tr1[i] == ',')
			{
			  tr1 = StringFunctions.ChangeCharacter(tr1, i, ':');
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
			  //(*testout) << "found rule" << endl;
			  netrule rule = new netrule();
			  //(*testout) << "fr1" << endl;
			  rule.LoadRule(ist);
			  //(*testout) << "fr2" << endl;
        
			  rules.Append(rule);
			  }
			  //(*testout) << "loop" << endl;
		  }
		  //(*testout) << "POS3" << endl;
        
		  ist = null;
		  //delete [] tr1;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public int ApplyRules(Array<Point2d> lpoints, Array<int> legalpoints, int maxlegalpoint, Array<INDEX_2> llines1, int maxlegalline, Array<Element2d> elements, Array<int> dellines, int tolerance, MeshingParameters mp)
		  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("meshing2::ApplyRules");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(ApplyRules_timer);
        
        
        
			double maxerr = 0.5 + 0.3 * tolerance;
			double minelerr = 2 + 0.5 * tolerance * tolerance;
        
			int noldlp = lpoints.Size();
			int noldll = llines1.Size();
        
        
			ArrayMem<int,100> pused = new ArrayMem<int,100>((uint)maxlegalpoint);
			ArrayMem<int,100> lused = new ArrayMem<int,100>((uint)maxlegalline);
			ArrayMem<int,100> pnearness = new ArrayMem<int,100>((uint)noldlp);
			ArrayMem<int,100> lnearness = new ArrayMem<int,100>(llines1.Size());
        
			ArrayMem<int, 20> pmap = new ArrayMem<int, 20>();
			ArrayMem<int, 20> pfixed = new ArrayMem<int, 20>();
			ArrayMem<int, 20> lmap = new ArrayMem<int, 20>();
        
			ArrayMem<Point2d,100> tempnewpoints = new ArrayMem<Point2d,100>();
			ArrayMem<INDEX_2,100> tempnewlines = new ArrayMem<INDEX_2,100>();
			ArrayMem<int,100> tempdellines = new ArrayMem<int,100>();
			ArrayMem<Element2d,100> tempelements = new ArrayMem<Element2d,100>();
        
        
			elements.SetSize(0);
			dellines.SetSize(0);
        
			testmode = debugparam.debugoutput;
        
		#if LOCDEBUG
			int loctestmode = testmode;
        
			if (loctestmode != 0)
			{
			(*testout) << "\n" << "\n" << "Check new environment" << "\n";
			(*testout) << "tolerance = " << tolerance << "\n";
			for (int i = 1; i <= lpoints.Size(); i++)
			{
			  (*testout) << "P" << i << " = " << lpoints.Get(i) << "\n";
			}
			(*testout) << "\n";
			for (int i = 1; i <= llines1.Size(); i++)
			{
			  (*testout) << "(" << llines1.Get(i).I1() << "-" << llines1.Get(i).I2() << ")" << "\n";
			}
			}
		#endif
        
			// check every rule
        
			int found = 0; // rule number
        
			pnearness = 1000;
        
			for (int j = 0; j < 2; j++)
			{
			  pnearness.Set(llines1[0][j], 0);
			}
        
        
        
        
			for (int cnt = 0; cnt < MAX_NEARNESS; cnt++)
			{
			bool ok = true;
			for (int i = 0; i < maxlegalline; i++)
			{
				INDEX_2 hline = llines1[i];
        
				int minn = min2(pnearness.Get(hline[0]), pnearness.Get(hline[1]));
        
				for (int j = 0; j < 2; j++)
				{
				  if (pnearness.Get(hline[j]) > minn + 1)
				  {
				  ok = false;
				  pnearness.Set(hline[j], minn + 1);
				  }
				}
			}
			if (!ok)
			{
				break;
			}
			}
        
        
			for (int i = 0; i < maxlegalline; i++)
			{
			  lnearness[i] = pnearness.Get(llines1[i][0]) + pnearness.Get(llines1[i][1]);
			}
        
        
			// resort lines after lnearness
			Array<INDEX_2> llines = new Array<INDEX_2>(llines1.Size());
			Array<int> sortlines = new Array<int>(llines1.Size());
			int[] lnearness_class = new int[MAX_NEARNESS];
        
			for (int j = 0; j < MAX_NEARNESS; j++)
			{
			  lnearness_class[j] = 0;
			}
			for (int i = 0; i < maxlegalline; i++)
			{
			  if (lnearness[i] < MAX_NEARNESS)
			  {
			lnearness_class[lnearness[i]]++;
			  }
			}
        
			int cumm = 0;
			for (int j = 0; j < MAX_NEARNESS; j++)
			{
			int hcnt = lnearness_class[j];
			lnearness_class[j] = cumm;
			cumm += hcnt;
			}
        
			for (int i = 0; i < maxlegalline; i++)
			{
			  if (lnearness[i] < MAX_NEARNESS)
			  {
			  llines[lnearness_class[lnearness[i]]] = llines1[i];
			  sortlines[lnearness_class[lnearness[i]]] = i + 1;
			  lnearness_class[lnearness[i]]++;
			  }
			  else
			  {
			  llines[cumm] = llines1[i];
			  sortlines[cumm] = i + 1;
			  cumm++;
			  }
			}
        
			for (int i = maxlegalline; i < llines1.Size(); i++)
			{
			llines[cumm] = llines1[i];
			sortlines[cumm] = i + 1;
			cumm++;
			}
        
			for (int i = 0; i < maxlegalline; i++)
			{
			  lnearness[i] = pnearness.Get(llines[i][0]) + pnearness.Get(llines[i][1]);
			}
        
        
        
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static bool firsttime = true;
			// static int timers[100];
			// static int timers2[100];
			// static int timers3[100];
			if (ApplyRules_firsttime)
			{
			/*
			for (int ri = 0; ri < rules.Size(); ri++)
			  timers[ri] = NgProfiler::CreateTimer (string("netrule ")+rules[ri]->Name());
			for (int ri = 0; ri < rules.Size(); ri++)
			  timers2[ri] = NgProfiler::CreateTimer (string("netrule,mapped ")+rules[ri]->Name());
			for (int ri = 0; ri < rules.Size(); ri++)
			  timers3[ri] = NgProfiler::CreateTimer (string("netrule,lines mapped ")+rules[ri]->Name());
			*/
			ApplyRules_firsttime = false;
			}
        
			lused = 0;
			pused = 0;
        
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1 = NgProfiler::CreateTimer("meshing2::ApplyRules 1");
			NgProfiler.RegionTimer reg1 = new NgProfiler.RegionTimer(ApplyRules_timer1);
        
        
			for (int ri = 1; ri <= rules.Size(); ri++)
			{
			// NgProfiler::RegionTimer reg(timers[ri-1]);
			netrule rule = rules.Get(ri);
        
		#if LOCDEBUG
			if (loctestmode != 0)
			{
			  (*testout) << "Rule " << rule.Name() << "\n";
			}
		#endif
        
			if (rule.GetQuality() > tolerance)
			{
				continue;
			}
        
			pmap.SetSize(rule.GetNP());
			lmap.SetSize(rule.GetNL());
        
			pmap = 0;
			lmap = 0;
        
			lused[0] = 1;
			lmap[0] = 1;
        
			for (int j = 0; j < 2; j++)
			{
				pmap.Elem(rule.GetLine 1[j]) = llines[0][j];
				pused.Elem(llines[0][j])++;
			}
        
        
        
			int nlok = 2;
        
        
			bool ok = false;
        
			while (nlok >= 2)
			{
        
				if (nlok <= rule.GetNOldL())
        
				{
				ok = false;
        
				int maxline = (rule.GetLNearness(nlok) < MAX_NEARNESS) ? lnearness_class[rule.GetLNearness(nlok)] : maxlegalline;
				// int maxline = maxlegalline;
        
				while (!ok && lmap.Get(nlok) < maxline)
				{
					lmap.Elem(nlok)++;
					int locli = lmap.Get(nlok);
        
					if (lnearness.Get(locli) > rule.GetLNearness(nlok))
					{
						continue;
					}
					if (lused.Get(locli))
					{
						continue;
					}
        
        
					ok = true;
        
					INDEX_2 loclin = llines.Get(locli);
					Vec2d linevec = lpoints.Get(loclin.I2()) - lpoints.Get(loclin.I1());
        
					if (rule.CalcLineError(nlok, linevec) > maxerr)
					{
					ok = false;
		#if LOCDEBUG
					if (loctestmode != 0)
					{
					  (*testout) << "not ok pos1" << "\n";
					}
		#endif
					continue;
					}
        
					for (int j = 0; j < 2; j++)
					{
					int refpi = rule.GetLine(nlok)[j];
        
					if (pmap.Get(refpi) != 0)
					{
						if (pmap.Get(refpi) != loclin[j])
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "not ok pos2" << "\n";
						}
		#endif
						break;
						}
					}
					else
					{
						if (rule.CalcPointDist(refpi, lpoints.Get(loclin[j])) > maxerr || !legalpoints.Get(loclin[j]) || pused.Get(loclin[j]))
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
							(*testout) << "nok pos3" << "\n";
							//if(rule->CalcPointDist (refpi, lpoints.Get(loclin[j])) > maxerr)
							//(*testout) << "r1" << endl;
							//if(!legalpoints.Get(loclin[j]))
							//(*testout) << "r2 legalpoints " << legalpoints << " loclin " << loclin << " j " << j << endl;
							//if(pused.Get(loclin[j]))
							//(*testout) << "r3" << endl;
						}
		#endif
						break;
						}
					}
					}
				}
        
				if (ok)
				{
					int locli = lmap.Get(nlok);
					INDEX_2 loclin = llines.Get(locli);
        
					lused.Elem(locli) = 1;
					for (int j = 0; j < 2; j++)
					{
					pmap.Set(rule.GetLine nlok[j], loclin[j]);
					pused.Elem(loclin[j])++;
					}
        
					nlok++;
				}
				else
				{
					lmap.Elem(nlok) = 0;
					nlok--;
        
					lused.Elem(lmap.Get(nlok)) = 0;
					for (int j = 0; j < 2; j++)
					{
					pused.Elem(llines.Get[] lmap.Get(nlok))--;
					if (!pused.Get(llines.Get[] lmap.Get(nlok)))
					{
					  pmap.Set(rule.GetLine nlok[j], 0);
					}
					}
				}
				}
        
				else
        
				{
				// NgProfiler::RegionTimer reg(timers3[ri-1]);
        
				// all lines are mapped !!
        
				// map also all points:
        
				int npok = 1;
				int incnpok = 1;
        
				pfixed.SetSize(pmap.Size());
				for (int i = 0; i < pmap.Size(); i++)
				{
				  pfixed[i] = (pmap[i] >= 1);
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
						ok = false;
        
						if (pmap.Get(npok))
						{
						  pused.Elem(pmap.Get(npok))--;
						}
        
						while (!ok && pmap.Get(npok) < maxlegalpoint)
						{
						ok = true;
        
						pmap.Elem(npok)++;
        
						if (pused.Get(pmap.Get(npok)))
						{
							ok = false;
						}
						else
						{
							if (rule.CalcPointDist(npok, lpoints.Get(pmap.Get(npok))) > maxerr || !legalpoints.Get(pmap.Get(npok)))
							{
        
							  ok = false;
							}
						}
						}
        
						if (ok)
						{
						pused.Elem(pmap.Get(npok))++;
						npok++;
						incnpok = 1;
						}
        
						else
        
						{
						pmap.Elem(npok) = 0;
						npok--;
						incnpok = 0;
						}
					}
					}
        
					else
        
					{
					// NgProfiler::RegionTimer reg(timers2[ri-1]);
        
					npok = rule.GetNOldP();
					incnpok = 0;
        
					if (ok)
					{
					  foundmap.Elem(ri)++;
					}
        
		#if LOCDEBUG
					if (loctestmode != 0)
					{
					  (*testout) << "lines and points mapped" << "\n";
					}
		#endif
        
					ok = true;
        
					// check orientations
        
					for (int i = 1; i <= rule.GetNOrientations(); i++)
					{
						if (CW(lpoints.Get(pmap.Get(rule.GetOrientation(i).i1)), lpoints.Get(pmap.Get(rule.GetOrientation(i).i2)), lpoints.Get(pmap.Get(rule.GetOrientation(i).i3))))
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "Orientation " << i << " not ok" << "\n";
						}
		#endif
						break;
						}
					}
        
        
					if (!ok)
					{
						continue;
					}
        
					Vector oldu = new Vector(2 * rule.GetNOldP());
        
					for (int i = 1; i <= rule.GetNOldP(); i++)
					{
						Vec2d ui = new Vec2d(rule.GetPoint(i), lpoints.Get(pmap.Get(i)));
						oldu(2 * i - 2) = ui.X();
						oldu(2 * i - 1) = ui.Y();
					}
        
					rule.SetFreeZoneTransformation(oldu, tolerance);
        
        
					if (!ok)
					{
						continue;
					}
					if (rule.ConvexFreeZone() == 0)
					{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "freezone not convex" << "\n";
						}
		#endif
						/*
						  static int cnt = 0;
						  cnt++;
						  if (cnt % 100 == 0)
						  {
						  cout << "freezone not convex, cnt = " << cnt << "; rule = " << rule->Name() << endl;
						  (*testout) << "freezone not convex, cnt = " << cnt << "; rule = " << rule->Name() << endl;
						  (*testout) << "tol = " << tolerance << endl;
						  (*testout) << "maxerr = " << maxerr << "; minerr = " << minelerr << endl;
						  (*testout) << "freezone = " << rule->GetTransFreeZone() << endl;
						  }
						*/
					}
        
					// check freezone:
					if (!ok)
					{
						continue;
					}
					for (int i = 1; i <= maxlegalpoint && ok; i++)
					{
						if (!pused.Get(i) && rule.IsInFreeZone(lpoints.Get(i)))
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "Point " << i << " in freezone" << "\n";
						}
		#endif
						break;
						}
					}
        
					if (!ok)
					{
						continue;
					}
					for (int i = maxlegalpoint + 1; i <= lpoints.Size(); i++)
					{
						if (rule.IsInFreeZone(lpoints.Get(i)))
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "Point " << i << " in freezone" << "\n";
						}
		#endif
						break;
						}
					}
        
        
					if (!ok)
					{
						continue;
					}
					for (int i = 1; i <= maxlegalline; i++)
					{
						if (!lused.Get(i) && rule.IsLineInFreeZone(lpoints.Get(llines.Get(i).I1()), lpoints.Get(llines.Get(i).I2())) != 0)
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "line " << llines.Get(i).I1() << "-" << llines.Get(i).I2() << " in freezone" << "\n";
						}
		#endif
						break;
						}
					}
        
					if (!ok)
					{
						continue;
					}
        
					for (int i = maxlegalline+1; i <= llines.Size(); i++)
					{
						if (rule.IsLineInFreeZone(lpoints.Get(llines.Get(i).I1()), lpoints.Get(llines.Get(i).I2())) != 0)
						{
						ok = false;
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "line " << llines.Get(i).I1() << "-" << llines.Get(i).I2() << " in freezone" << "\n";
						}
		#endif
						break;
						}
					}
        
        
					/*
					// check orientations
		
					for (i = 1; i <= rule->GetNOrientations() && ok; i++)
					{
					if (CW (lpoints.Get(pmap.Get(rule->GetOrientation(i).i1)),
					lpoints.Get(pmap.Get(rule->GetOrientation(i).i2)),
					lpoints.Get(pmap.Get(rule->GetOrientation(i).i3))) )
					{
					ok = 0;
					if (loctestmode)
					(*testout) << "Orientation " << i << " not ok" << endl;
					}
					}
					*/
        
        
					if (!ok)
					{
						continue;
					}
        
		#if LOCDEBUG
					if (loctestmode != 0)
					{
					  (*testout) << "rule ok" << "\n";
					}
		#endif
        
					// Setze neue Punkte:
					if (rule.GetNOldP() < rule.GetNP())
					{
						Vector newu = new Vector(rule.GetOldUToNewU.functorMethod().Height());
						rule.GetOldUToNewU.functorMethod().Mult(oldu, newu);
        
						int oldnp = rule.GetNOldP();
						for (int i = oldnp + 1; i <= rule.GetNP(); i++)
						{
						Point2d np = rule.GetPoint(i);
						np.X() += newu(2 * (i - oldnp) - 2);
						np.Y() += newu(2 * (i - oldnp) - 1);
        
										lpoints.Append(np);
						pmap.Elem(i) = lpoints.Size();
						}
					}
        
					// Setze neue Linien:
        
					for (int i = rule.GetNOldL() + 1; i <= rule.GetNL(); i++)
					{
						llines.Append(new INDEX_2(pmap.Get(rule.GetLine i[0]), pmap.Get(rule.GetLine i[1])));
					}
        
        
					// delete old lines:
					for (int i = 1; i <= rule.GetNDelL(); i++)
					{
					  dellines.Append(sortlines.Elem(lmap.Get(rule.GetDelLine(i))));
					}
					// dellines.Append (lmap.Get(rule->GetDelLine(i))));
        
					// dellines.Append (lmap.Elem(rule->GetDelLines()));
					// lmap[rule->GetDelLines()];
        
        
					// insert new elements:
        
					for (int i = 1; i <= rule.GetNE(); i++)
					{
						elements.Append(rule.GetElement(i));
						for (int j = 1; j <= elements.Get(i).GetNP(); j++)
						{
						  elements.Elem(i).PNum(j) = pmap.Get(elements.Get(i).PNum(j));
						}
					}
        
        
					double elerr = 0;
					for (int i = 1; i <= elements.Size(); i++)
					{
						double hf;
						if (mp.quad == 0)
						{
						  hf = CalcElementBadness(lpoints, elements.Get(i));
						}
						else
						{
						  hf = elements.Get(i).CalcJacobianBadness(lpoints) * 5;
						}
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						  (*testout) << "r " << rule.Name() << "bad = " << hf << "\n";
						}
		#endif
						if (hf > elerr)
						{
							elerr = hf;
						}
					}
        
		#if LOCDEBUG
					if (loctestmode != 0)
					{
					  (*testout) << "error = " << elerr;
					}
		#endif
        
					canuse.Elem(ri)++;
        
					if (elerr < 0.99 * minelerr)
					{
		#if LOCDEBUG
						if (loctestmode != 0)
						{
						(*testout) << "rule = " << rule.Name() << "\n";
						(*testout) << "class = " << tolerance << "\n";
						(*testout) << "lpoints: " << "\n";
						for (int i = 1; i <= lpoints.Size(); i++)
						{
						  (*testout) << lpoints.Get(i) << "\n";
						}
						(*testout) << "llines: " << "\n";
						for (int i = 1; i <= llines.Size(); i++)
						{
						  (*testout) << llines.Get(i).I1() << " " << llines.Get(i).I2() << "\n";
						}
        
						(*testout) << "Freezone: ";
						for (int i = 1; i <= rule.GetTransFreeZone().Size(); i++)
						{
						  (*testout) << rule.GetTransFreeZone().Get(i) << "\n";
						}
						}
		#endif
        
						minelerr = elerr;
						found = ri;
        
						tempnewpoints = lpoints.Range(noldlp, lpoints.Size());
						tempnewlines = llines.Range(noldll, llines.Size());
						tempdellines = dellines;
						tempelements = elements;
					}
        
					lpoints.SetSize(noldlp);
					llines.SetSize(noldll);
					dellines.SetSize(0);
					elements.SetSize(0);
					ok = false;
					}
				}
        
				nlok = rule.GetNOldL();
        
				lused.Set(lmap.Get(nlok), 0);
        
				for (int j = 1; j <= 2; j++)
				{
					int refpi = rule.GetPointNr(nlok, j);
					pused.Elem(pmap.Get(refpi))--;
        
					if (pused.Get(pmap.Get(refpi)) == 0)
					{
					  pmap.Set(refpi, 0);
					}
				}
				}
			}
			}
        
        
			if (found != 0)
			{
			lpoints.Append(tempnewpoints);
			llines1.Append(tempnewlines);
			dellines.Append(tempdellines);
			elements.Append(tempelements);
			}
        
        
			return found;
		  }
	}
}