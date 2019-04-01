using System;
using System.Collections.Generic;

public static class GlobalMembers
{
	public static ostream operator << (ostream s, MiniElement2d el)
	{
	  s << "np = " << el.GetNP();
	  for (int j = 0; j < el.GetNP(); j++)
	  {
		s << " " << el[j];
	  }
	  return s;
	}



//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//DLL_HEADER void BisectTetsCopyMesh(Mesh UnnamedParameter, NetgenGeometry UnnamedParameter2, BisectionOptions opt);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//DLL_HEADER void ZRefinement(Mesh UnnamedParameter, NetgenGeometry UnnamedParameter2, ZRefinementOptions opt);
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_EXPORT __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_IMPORT __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_EXPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_IMPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max) { if ((value)<(min) || (value)>=(max)) throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max)); }
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_START(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_END(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_TRACER(n)


	///
	public static void InsertVirtualBoundaryLayer(Mesh mesh)
	{
	   Console.Write("Insert virt. b.l.");
	   Console.Write("\n");

	   int surfid;

	   Console.Write("Boundary Nr:");
	   surfid = int.Parse(ConsoleInput.ReadToWhiteSpace(true));

	   int i;
	   int np = mesh.GetNP();

	   Console.Write("Old NP: ");
	   Console.Write(mesh.GetNP());
	   Console.Write("\n");
	   Console.Write("Trigs: ");
	   Console.Write(mesh.GetNSE());
	   Console.Write("\n");

	   BitArray bndnodes = new BitArray(np);
	   Array<int> mapto = new Array<int>(np);

	   bndnodes.Clear();
	   for (i = 1; i <= mesh.GetNSeg(); i++)
	   {
		  int snr = mesh.LineSegment(i).edgenr;
		  Console.Write("snr = ");
		  Console.Write(snr);
		  Console.Write("\n");
		  if (snr == surfid)
		  {
			 bndnodes.Set(mesh.LineSegment i[0]);
			 bndnodes.Set(mesh.LineSegment i[1]);
		  }
	   }
	   for (i = 1; i <= mesh.GetNSeg(); i++)
	   {
		  int snr = mesh.LineSegment(i).edgenr;
		  if (snr != surfid)
		  {
			 bndnodes.Clear(mesh.LineSegment i[0]);
			 bndnodes.Clear(mesh.LineSegment i[1]);
		  }
	   }

	   for (i = 1; i <= np; i++)
	   {
		   if (bndnodes.Test(i))
		   {
			 mapto.Elem(i) = mesh.AddPoint(new mesh.Point(i));
		   }
		   else
		   {
			 mapto.Elem(i) = 0;
		   }
	   }

	   for (i = 1; i <= mesh.GetNSE(); i++)
	   {
		  Element2d el = mesh.SurfaceElement(i);
		  for (int j = 1; j <= el.GetNP(); j++)
		  {
			 if (mapto.Get(el.PNum(j)))
			 {
				el.PNum(j) = mapto.Get(el.PNum(j));
			 }
		  }
	   }


	   int nq = 0;
	   for (i = 1; i <= mesh.GetNSeg(); i++)
	   {
		  int snr = mesh.LineSegment(i).edgenr;
		  if (snr == surfid)
		  {
			 int p1 = mesh.LineSegment(i)[0];
			 int p2 = mesh.LineSegment(i)[1];
			 int p3 = mapto.Get(p1);
			 if (p3 == 0)
			 {
				 p3 = p1;
			 }
			 int p4 = mapto.Get(p2);
			 if (p4 == 0)
			 {
				 p4 = p2;
			 }

			 Element2d el = new Element2d(QUAD);
			 el.PNum(1) = p1;
			 el.PNum(2) = p2;
			 el.PNum(3) = p3;
			 el.PNum(4) = p4;
			 el.SetIndex(2);
			 mesh.AddSurfaceElement(el);
			 nq++;
		  }
	   }

	   Console.Write("New NP: ");
	   Console.Write(mesh.GetNP());
	   Console.Write("\n");
	   Console.Write("Quads: ");
	   Console.Write(nq);
	   Console.Write("\n");
	}

/*
    Philippose Rajan - 11 June 2009
    
    Added an initial experimental function for 
    generating prismatic boundary layers on 
    a given set of surfaces.
    
    The number of layers, height of the first layer 
    and the growth / shrink factor can be specified 
    by the user

    Currently, the layer height is calculated using:
    height = h_first_layer * (growth_factor^(num_layers - 1))
*/

	public static void GenerateBoundaryLayer(Mesh mesh, BoundaryLayerParameters blp)
	{
	   ofstream dbg = new ofstream("BndLayerDebug.log");

	   // Angle between a surface element and a growth-vector below which
	   // a prism is project onto that surface as a quad
	   // (in degrees)
	   double angleThreshold = 5.0;


	   Array<int> surfid = new Array<int>(blp.surfid);
	   int prismlayers = blp.prismlayers;
	   double hfirst = blp.hfirst;
	   double growthfactor = blp.growthfactor;
	   Array<double> heights = new Array<double>(blp.heights);

	   bool grow_edges = false; // grow layer at edges


	   // Monitor and print out the number of prism and quad elements
	   // added to the mesh
	   int numprisms = 0;
	   int numquads = 0;


	   Console.Write("Old NP: ");
	   Console.Write(mesh.GetNP());
	   Console.Write("\n");
	   Console.Write("Old NSE: ");
	   Console.Write(mesh.GetNSE());
	   Console.Write("\n");

	   for (int layer = prismlayers; layer >= 1; layer--)
	   {
		   Console.Write("Generating layer: ");
		   Console.Write(layer);
		   Console.Write("\n");

		   MeshTopology meshtopo = mesh.GetTopology();
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
		   const_cast<MeshTopology &> (meshtopo).SetBuildEdges(true);
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
		   const_cast<MeshTopology &> (meshtopo).SetBuildFaces(true);
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
		   const_cast<MeshTopology &> (meshtopo).Update();

		   double layerht = hfirst;

		   if (heights.Size() > 0)
		   {
			   layerht = heights[layer - 1];
		   }
		   else
		   {
			   if (growthfactor == 1)
			   {
				   layerht = layer * hfirst;
			   }
			   else
			   {
				   layerht = hfirst * (ngsimd.GlobalMembers.pow(growthfactor, (layer + 1)) - 1) / (growthfactor - 1);
			   }
		   }

		  Console.Write("Layer Height = ");
		  Console.Write(layerht);
		  Console.Write("\n");

		  // Need to store the old number of points and
		  // surface elements because there are new points and
		  // surface elements being added during the process
		  int np = mesh.GetNP();
		  int nse = mesh.GetNSE();
		  int ne = mesh.GetNE();

		  // Safety measure to ensure no issues with mesh
		  // consistency
		  int nseg = mesh.GetNSeg();

		  // Indicate which points need to be remapped
		  BitArray bndnodes = new BitArray(np + 1); // big enough for 1-based array

		  // Map of the old points to the new points
		  Array<PointIndex, PointIndex.BASE> mapto = new Array<PointIndex, PointIndex.BASE>(np);

		  // Growth vectors for the prismatic layer based on
		  // the effective surface normal at a given point
		  Array<Vec3d, PointIndex.BASE> growthvectors = new Array<Vec3d, PointIndex.BASE>(np);

		  // Bit array to identify all the points belonging
		  // to the surface of interest
		  bndnodes.Clear();

		  // Run through all the surface elements and mark the points
		  // belonging to those where a boundary layer has to be created.
		  // In addition, also calculate the effective surface normal
		  // vectors at each of those points to determine the mesh motion
		  // direction
		  Console.Write("Marking points for remapping....");
		  Console.Write("\n");

		  for (SurfaceElementIndex si = 0; si < nse; si++)
		  {
			if (surfid.Contains(mesh[si].GetIndex()))
			{
				Element2d sel = mesh[si];
				for (int j = 0; j < sel.GetNP(); j++)
				{
					// Set the bitarray to indicate that the
					// point is part of the required set
					bndnodes.Set(sel[j]);
					Vec3d surfacenormal = new Vec3d();

					// Calculate the surface normal at the current point
					// with respect to the current surface element
					netgen.GlobalMembers.GetSurfaceNormal(mesh,sel,j + 1,surfacenormal);

					// Add the surface normal to the already existent one
					// (This gives the effective normal direction at corners
					//  and curved areas)
					growthvectors[sel[j]] += surfacenormal;
				}
			}
		  }

		  if (!grow_edges)
		  {
			for (SegmentIndex sei = 0; sei <= nseg; sei++)
			{
				bndnodes.Clear(mesh[sei][0]);
				bndnodes.Clear(mesh[sei][1]);
			}
		  }

		  // Add additional points into the mesh structure in order to
		  // clone the surface elements.
		  // Also invert the growth vectors so that they point inwards,
		  // and normalize them
		  Console.Write("Cloning points and calculating growth vectors....");
		  Console.Write("\n");

		  for (PointIndex pi = 1; pi <= np; pi++)
		  {
			  if (bndnodes.Test(pi))
			  {
				  mapto[pi] = mesh.AddPoint(mesh[pi]);

				  growthvectors[pi].Normalize();
				  growthvectors[pi] *= -1.0;
			  }
			  else
			  {
				  mapto[pi] = 0;
				  growthvectors[pi] = new Vec3d(0, 0, 0);
			  }
		  }


		  // Add quad surface elements at edges for surfaces which
		  // don't have boundary layers

		  // Bit array to keep track of segments already processed
		  BitArray segsel = new BitArray(nseg);

		  // Set them all to "1" to initially activate all segments
		  segsel.Set();

		  Console.Write("Adding 2D Quad elements on required surfaces....");
		  Console.Write("\n");

		  if (grow_edges)
		  {
		  for (SegmentIndex sei = 0; sei <= nseg; sei++)
		  {
			  PointIndex seg_p1 = mesh[sei][0];
			  PointIndex seg_p2 = mesh[sei][1];

			  // Only go in if the segment is still active, and if both its
			  // surface index is part of the "hit-list"
			  if (segsel.Test(sei) && surfid.Contains(mesh[sei].si))
			  {
				  // clear the bit to indicate that this segment has been processed
				  segsel.Clear(sei);

				  // Find matching segment pair on other surface
				  for (SegmentIndex sej = 0; sej < nseg; sej++)
				  {
					  PointIndex segpair_p1 = mesh[sej][1];
					  PointIndex segpair_p2 = mesh[sej][0];

					  // Find the segment pair on the neighbouring surface element
					  // Identified by: seg1[0] = seg_pair[1] and seg1[1] = seg_pair[0]
					  if (segsel.Test(sej) && ((segpair_p1 == seg_p1) && (segpair_p2 == seg_p2)))
					  {
						  // clear bit to indicate that processing of this segment is done
						  segsel.Clear(sej);

						  // Only worry about those surfaces which are not in the
						  // boundary layer list
						  if (!surfid.Contains(mesh[sej].si))
						  {
							  SurfaceElementIndex pnt_commelem = 0;
							  Array<SurfaceElementIndex> pnt1_elems = new Array<SurfaceElementIndex>();
							  Array<SurfaceElementIndex> pnt2_elems = new Array<SurfaceElementIndex>();


							  meshtopo.GetVertexSurfaceElements(segpair_p1,pnt1_elems);
							  meshtopo.GetVertexSurfaceElements(segpair_p2,pnt2_elems);

							  for (int k = 0; k < pnt1_elems.Size(); k++)
							  {
								  Element2d pnt1_sel = mesh.SurfaceElement(pnt1_elems[k]);
								  for (int l = 0; l < pnt2_elems.Size(); l++)
								  {
									  Element2d pnt2_sel = mesh.SurfaceElement(pnt2_elems[l]);
									  if ((pnt1_sel.GetIndex() == mesh[sej].si) && (pnt2_sel.GetIndex() == mesh[sej].si) && (pnt1_elems[k] == pnt2_elems[l]))
									  {
										  pnt_commelem = pnt1_elems[k];
									  }
								  }
							  }

							  /*
							    int pnum_commelem = 0;
							    for(int k = 1; k <= mesh.SurfaceElement(pnt_commelem).GetNP(); k++)
							    {
							    if((mesh.SurfaceElement(pnt_commelem).PNum(k) != segpair_p1)
							    && (mesh.SurfaceElement(pnt_commelem).PNum(k) != segpair_p2))
							    {
							    pnum_commelem = mesh.SurfaceElement(pnt_commelem).PNum(k);
							    }
							    }
							  */

							  Vec3d surfelem_vect = new Vec3d();
							  Vec3d surfelem_vect1 = new Vec3d();

							  Element2d commsel = mesh.SurfaceElement(pnt_commelem);

							  dbg << "NP= " << commsel.GetNP() << " : ";

							  for (int k = 1; k <= commsel.GetNP(); k++)
							  {
								  netgen.GlobalMembers.GetSurfaceNormal(mesh,commsel,k,surfelem_vect1);
								  surfelem_vect += surfelem_vect1;
							  }

							  surfelem_vect.Normalize();

							  double surfangle = Angle(growthvectors.Elem(segpair_p1),surfelem_vect);

							  dbg << "V1= " << surfelem_vect1 << " : V2= " << surfelem_vect1 << " : V= " << surfelem_vect << " : GV= " << growthvectors.Elem(segpair_p1) << " : Angle= " << surfangle * 180 / 3.141592;


							  // remap the segments to the new points
							  mesh[sei][0] = mapto[seg_p1];
							  mesh[sei][1] = mapto[seg_p2];
							  mesh[sej][1] = mapto[seg_p1];
							  mesh[sej][0] = mapto[seg_p2];

							  if ((surfangle < (90 + angleThreshold) * 3.141592 / 180.0) && (surfangle > (90 - angleThreshold) * 3.141592 / 180.0))
							  {
								  dbg << " : quad\n";
								  // Since the surface is lower than the threshold, change the effective
								  // prism growth vector to match with the surface vector, so that
								  // the Quad which is created lies on the original surface
								  //growthvectors.Elem(segpair_p1) = surfelem_vect;

								  // Add a quad element to account for the prism volume
								  // element which is going to be added
								  Element2d sel = new Element2d(QUAD);
								  sel.PNum(4) = mapto[seg_p1];
								  sel.PNum(3) = mapto[seg_p2];
								  sel.PNum(2) = segpair_p2;
								  sel.PNum(1) = segpair_p1;
								  sel.SetIndex(mesh[sej].si);
								  mesh.AddSurfaceElement(sel);
								  numquads++;
							  }
							  else
							  {
								  dbg << "\n";
								  for (int k = 0; k < pnt1_elems.Size(); k++)
								  {
									  Element2d pnt_sel = mesh.SurfaceElement(pnt1_elems[k]);
									  if (pnt_sel.GetIndex() == mesh[sej].si)
									  {
										  for (int l = 0; l < pnt_sel.GetNP(); l++)
										  {
											  if (pnt_sel[l] == segpair_p1)
											  {
												pnt_sel[l] = mapto[seg_p1];
											  }
											  else if (pnt_sel[l] == segpair_p2)
											  {
												pnt_sel[l] = mapto[seg_p2];
											  }
										  }
									  }
								  }

								  for (int k = 0; k < pnt2_elems.Size(); k++)
								  {
									  Element2d pnt_sel = mesh.SurfaceElement(pnt2_elems[k]);
									  if (pnt_sel.GetIndex() == mesh[sej].si)
									  {
										  for (int l = 0; l < pnt_sel.GetNP(); l++)
										  {
											  if (pnt_sel[l] == segpair_p1)
											  {
												pnt_sel[l] = mapto.Get(seg_p1);
											  }
											  else if (pnt_sel[l] == segpair_p2)
											  {
												pnt_sel[l] = mapto.Get(seg_p2);
											  }
										  }
									  }
								  }
							  }
							  // }
						  }
						  else
						  {
							  // If the code comes here, it indicates that we are at
							  // a line segment pair which is at the intersection
							  // of two surfaces, both of which have to grow boundary
							  // layers.... here too, remapping the segments to the
							  // new points is required
							  mesh[sei][0] = mapto.Get(seg_p1);
							  mesh[sei][1] = mapto.Get(seg_p2);
							  mesh[sej][1] = mapto.Get(seg_p1);
							  mesh[sej][0] = mapto.Get(seg_p2);
						  }
					  }
				  }
			  }
		  }
		  }

		  // Add prismatic cells at the boundaries
		  Console.Write("Generating prism boundary layer volume elements....");
		  Console.Write("\n");

		  for (SurfaceElementIndex si = 0; si < nse; si++)
		  {
			  Element2d sel = mesh.SurfaceElement(si);
			  if (surfid.Contains(sel.GetIndex()))
			  {
				  /*
				  Element el(PRISM);
				  for (int j = 0; j < sel.GetNP(); j++)
				    {
				      // Check (Doublecheck) if the corresponding point has a
				      // copy available for remapping
				      if (mapto.Get(sel[j]))
				        {
				          // Define the points of the newly added Prism cell
				          el[j+3] = mapto[sel[j]];
				          el[j] = sel[j];
				        }
				      else
				        {
				          el[j+3] = sel[j];
				          el[j] = sel[j];
				        }
				    }
				  
				  el.SetIndex(1);
				  el.Invert();
				  mesh.AddVolumeElement(el);
				  numprisms++;
				  */
				  // cout << "add element: " << endl;
				  int classify = 0;
				  for (int j = 0; j < 3; j++)
				  {
					if (mapto[sel[j]])
					{
					  classify += (1 << j);
					}
				  }

				  // cout << "classify = " << classify << endl;

				  ELEMENT_TYPE[] types = {PRISM, TET, TET, PYRAMID, TET, PYRAMID, PYRAMID, PRISM};
				  int[] nums = {sel[0], sel[1], sel[2], mapto[sel[0]], mapto[sel[1]], mapto[sel[2]]};
				  int[][] vertices =
				  {
					  new int[] {0, 1, 2, 0, 1, 2},
					  new int[] {0, 2, 1, 3, 0, 0},
					  new int[] {0, 2, 1, 4, 0, 0},
					  new int[] {0, 1, 4, 3, 2, 0},
					  new int[] {0, 2, 1, 5, 0, 0},
					  new int[] {2, 0, 3, 5, 1, 0},
					  new int[] {1, 2, 5, 4, 0, 0},
					  new int[] {0, 2, 1, 3, 5, 4}
				  };

				  Element el = new Element(types[classify]);
				  for (int i = 0; i < 6; i++)
				  {
					el[i] = nums[vertices[classify][i]];
				  }
					if (blp.new_matnrs.Size() > 0)
					{
					   el.SetIndex(blp.new_matnrs[layer - 1]);
					}
					else
					{
					   el.SetIndex(blp.new_matnr);
					}
				  // cout << "el = " << el << endl;
				  if (classify != 0)
				  {
					mesh.AddVolumeElement(el);
				  }
			  }
		  }

		  // Finally switch the point indices of the surface elements
		  // to the newly added ones
		  Console.Write("Transferring boundary layer surface elements to new vertex references....");
		  Console.Write("\n");

		  for (int i = 1; i <= nse; i++)
		  {
			  Element2d sel = mesh.SurfaceElement(i);
			  if (surfid.Contains(sel.GetIndex()))
			  {
				 for (int j = 1; j <= sel.GetNP(); j++)
				 {
					 // Check (Doublecheck) if the corresponding point has a
					 // copy available for remapping
					 if (mapto.Get(sel.PNum(j)))
					 {
						 // Map the surface elements to the new points
						 sel.PNum(j) = mapto.Get(sel.PNum(j));
					 }
				 }
			  }
		  }
		  for (int i = 1; i <= ne; i++)
		  {
			  Element el = mesh.VolumeElement(i);
			  if (el.GetIndex() != blp.bulk_matnr)
			  {
				 for (int j = 1; j <= el.GetNP(); j++)
				 {
					 // Check (Doublecheck) if the corresponding point has a
					 // copy available for remapping
					 if (mapto.Get(el.PNum(j)))
					 {
						 // Map the surface elements to the new points
						 el.PNum(j) = mapto.Get(el.PNum(j));
					 }
				 }
			  }
		  }




		  // Lock all the prism points so that the rest of the mesh can be
		  // optimised without invalidating the entire mesh
		  for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
		  {
			if (bndnodes.Test(pi))
			{
				mesh.AddLockedPoint(pi);
			}
		  }

		  // Now, actually pull back the old surface points to create
		  // the actual boundary layers
		  Console.Write("Moving and optimising boundary layer points....");
		  Console.Write("\n");

		  for (int i = 1; i <= np; i++)
		  {
			 Array<ElementIndex> vertelems = new Array<ElementIndex>();

			 if (bndnodes.Test(i))
			 {
				 MeshPoint pointtomove = new MeshPoint();

				 pointtomove = new mesh.Point(i);

				 if (layer == prismlayers)
				 {
					 mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors.Elem(i));

					 meshtopo.GetVertexElements(i,vertelems);

					 for (int j = 1; j <= vertelems.Size(); j++)
					 {
						 // double sfact = 0.9;
						 Element volel = mesh.VolumeElement(vertelems.Elem(j));
						 if (((volel.GetType() == TET) || (volel.GetType() == TET10)) && (!volel.IsDeleted()))
						 {
							 //while((volel.Volume(mesh.Points()) <= 0.0) && (sfact >= 0.0))
							 //{
							 //   mesh.Point(i).SetPoint(pointtomove + (sfact * layerht * growthvectors.Elem(i)));
							 //   mesh.ImproveMesh();

							 //   // Try to move the point back by one step but
							 //   // if the volume drops to below zero, double back
							 //   mesh.Point(i).SetPoint(pointtomove + ((sfact + 0.1) * layerht * growthvectors.Elem(i)));
							 //   if(volel.Volume(mesh.Points()) <= 0.0)
							 //   {
							 //      mesh.Point(i).SetPoint(pointtomove + (sfact * layerht * growthvectors.Elem(i)));
							 //   }
							 //   sfact -= 0.1;
							 //}
							 volel.Delete();
						 }
					 }
				 }
				 else
				 {
					 mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors.Elem(i));
				 }
			 }
		  }
	  mesh.Compress();
	   }

	   // Optimise the tet part of the volume mesh after all the modifications
	   // to the system are completed
	   //OptimizeVolume(mparam,mesh);

	   Console.Write("New NP: ");
	   Console.Write(mesh.GetNP());
	   Console.Write("\n");
	   Console.Write("Num of Quads: ");
	   Console.Write(numquads);
	   Console.Write("\n");
	   Console.Write("Num of Prisms: ");
	   Console.Write(numprisms);
	   Console.Write("\n");
	   Console.Write("Boundary Layer Generation....Done!");
	   Console.Write("\n");

	   dbg.close();
	}

/* ***************************** HPRefinement ********************************** */




	public static void HPRefinement(Mesh mesh, Refinement @ref, int levels, double fac1 = 0.125, bool setorders = true, bool reflevels = false)
	{
	  PrintMessage(1, "HP Refinement called, levels = ", levels);


	  // NgLock mem_lock (mem_mutex,1);

	  mesh.coarsemesh = new Mesh();
	  *mesh.coarsemesh = mesh;

	  // #ifdef CURVEDELEMS_NEW
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
	  const_cast<CurvedElements&> (mesh.coarsemesh.GetCurvedElements()).BuildCurvedElements(@ref, mesh.GetCurvedElements().GetOrder());
	  // #endif


	  mesh.hpelements = null;
	  mesh.hpelements = new Array<HPRefElement>();

	  Array<HPRefElement> hpelements = *mesh.hpelements;

	  netgen.GlobalMembers.InitHPElements(mesh,hpelements);

	  Array<int> nplevel = new Array<int>();
	  nplevel.Append(mesh.GetNP());

	  int act_ref = 1;
	  bool sing = netgen.GlobalMembers.ClassifyHPElements(mesh,hpelements, act_ref, levels);

	  sing = true; // iterate at least once
	  while (sing)
	  {
	  Console.Write(" Start new hp-refinement: step ");
	  Console.Write(act_ref);
	  Console.Write("\n");

	  netgen.GlobalMembers.DoRefinement(mesh, hpelements, @ref, fac1);
	  netgen.GlobalMembers.DoRefineDummies(mesh, hpelements, @ref);

	  nplevel.Append(mesh.GetNP());
	  netgen.GlobalMembers.CalcStatistics(hpelements);

	  netgen.GlobalMembers.SubdivideDegeneratedHexes(mesh, hpelements,fac1);

		  netgen.GlobalMembers.ReorderPoints(mesh, hpelements);

	  mesh.ClearSegments();
	  mesh.ClearSurfaceElements();
		mesh.ClearVolumeElements();

	  for (int i = 0; i < hpelements.Size(); i++)
	  {
		  HPRefElement hpel = hpelements[i];
		  if (netgen.GlobalMembers.Get_HPRef_Struct(hpel.type))
		  {
			switch (netgen.GlobalMembers.Get_HPRef_Struct(hpel.type).geom)
			{
		  case HP_SEGM:
		  {
			  Segment seg = new Segment();
			  seg[0] = hpel.pnums[0];
			  seg[1] = hpel.pnums[1];
			  // NOTE: only for less than 10000 elements (HACK) !!!
			  seg.edgenr = hpel.index % 10000;
			  seg.si = hpel.index / 10000;

					  /*
					  seg.epgeominfo[0].dist = hpel.param[0][0]; // he: war hpel.param[0][0]
					  seg.epgeominfo[1].dist = hpel.param[1][0]; // he: war hpel.param[1][0]
					  */

					  Segment coarseseg = mesh.coarsemesh.LineSegment(hpel.coarse_elnr + 1);
					  double d1 = coarseseg.epgeominfo[0].dist;
					  double d2 = coarseseg.epgeominfo[1].dist;

					  // seg.epgeominfo[0].dist = hpel.param[0][0]; // he: war hpel.param[0][0]
					  // seg.epgeominfo[1].dist = hpel.param[1][0]; // he: war hpel.param[1][0]

					  seg.epgeominfo[0].dist = d1 + hpel.param[0][0] * (d2 - d1); // JS, June 08
					  seg.epgeominfo[1].dist = d1 + hpel.param[1][0] * (d2 - d1);


			  seg.epgeominfo[0].edgenr = seg.edgenr;
			  seg.epgeominfo[1].edgenr = seg.edgenr;
					  seg.domin = hpel.domin;
					  seg.domout = hpel.domout; // he: needed for segments!
			  seg.hp_elnr = i;
			  seg.singedge_left = hpel.singedge_left;
			  seg.singedge_right = hpel.singedge_right;
			  mesh.AddSegment(seg);
			  break;
		  }

		  case HP_TRIG:
		  case HP_QUAD:
		  {
			  Element2d el = new Element2d(hpel.np);
			  for (int j = 0;j < hpel.np;j++)
			  {
				el.PNum(j + 1) = hpel.pnums[j];
			  }
			  el.hp_elnr = i;
			  el.SetIndex(hpel.index);
			  if (setorders)
			  {
				el.SetOrder(act_ref + 1,act_ref + 1,0);
			  }
			  mesh.AddSurfaceElement(el);
			  break;
		  }
		  case HP_HEX:
		  case HP_TET:
		  case HP_PRISM:
		  case HP_PYRAMID:
		  {
			  Element el = new Element(hpel.np);
			  for (int j = 0;j < hpel.np;j++)
			  {
				el.PNum(j + 1) = hpel.pnums[j];
			  }
			  el.SetIndex(hpel.index);
			  el.hp_elnr = i;
			  if (setorders)
			  {
				el.SetOrder(act_ref + 1,act_ref + 1,act_ref + 1);
			  }
			  mesh.AddVolumeElement(el);
			  break;
		  }

		  default:
			PrintSysError("hpref, backconversion failed for element ", (int)(netgen.GlobalMembers.Get_HPRef_Struct(hpel.type).geom));
		break;
			}
		  }
	  }
	  Console.Write(" Start with Update Topology ");
	  Console.Write("\n");
	  mesh.UpdateTopology();
	  Console.Write(" Mesh Update Topology done ");
	  Console.Write("\n");

	  act_ref++;

	  sing = netgen.GlobalMembers.ClassifyHPElements(mesh,hpelements, act_ref, levels);
	  }

	  Console.Write(" HP-Refinement done with ");
	  Console.Write(--act_ref);
	  Console.Write(" refinement steps.");
	  Console.Write("\n");

	  if (act_ref >= 1)
	  {
	  for (ElementIndex i = 0;i < mesh.GetNE(); i++)
	  {
		  // Element el = mesh[i] ;
		  HPRefElement hpel = hpelements[mesh[i].hp_elnr];
		  ELEMENT_EDGE edges = MeshTopology.GetEdges1(mesh[i].GetType());
		  double[] dist = {0, 0, 0};
		  int[] ord_dir = {0, 0, 0};
		  int[] edge_dir = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		  int ned = 4;

		  switch (mesh[i].GetType())
		  {
			case TET:
		  /* cout << " TET " ;
		  for(int k=0;k<4;k++) cout << el[k] << "\t" ;
		  cout << endl; */
		  break;
			case PRISM:
		  /* cout << " PRISM " ;
		  for(int k=0;k<6;k++) cout << el[k] << "\t" ;
		  cout << endl;  */
		  for (int l = 6;l < 9;l++)
		  {
			  edge_dir[l] = 2;
		  }
		  ord_dir[2] = 2;
		  ned = 9;
		  break;
			case HEX:
		  /* cout << " HEX " ;
		  for(int k=0;k<8;k++) cout << el[k] << "\t" ;
		  cout << endl; */
		  for (int l = 8;l < 12; l++)
		  {
			  edge_dir[l] = 2;
		  }
		  edge_dir[2] = edge_dir[3] = edge_dir[6] = edge_dir[7] = 1;
		  ord_dir[1] = 1;
		  ord_dir[2] = 2;
		  ned = 12;
		  break;
			case PYRAMID:
		  /*	cout << " PYRAMID " ;
		  for(int k=0;k<5;k++) cout << el[k] << "\t" ;
		  cout << endl; */
		  for (int l = 4;l < 8;l++)
		  {
			  edge_dir[l] = 2;
		  }
		  edge_dir[2] = edge_dir[3] = 1;
		  ord_dir[1] = 1;
		  ord_dir[2] = 2;
		  ned = 8;
		  break;


				default:
				  cerr << "HPRefElement: illegal elementtype (2) " << mesh[i].GetType() << "\n";
				  throw new Exception("HPRefElement: illegal elementtype (2)");

		  }

		  for (int j = 0;j < ned;j++)
		  {

		  Vec < 3> v(hpel.param[edges[j][0] - 1][0] - hpel.param[edges[j][1] - 1][0], hpel.param[edges[j][0] - 1][1] - hpel.param[edges[j][1] - 1][1], hpel.param[edges[j][0] - 1][2] - hpel.param[edges[j][1] - 1][2]);
		  dist[edge_dir[j]] = Math.Max(v.Length(),dist[edge_dir[j]]);
		  }

		  int[] refi = new int[3];
		  for (int j = 0;j < 3;j++)
		  {
			refi[j] = (int)Math.Max((double)ngsimd.GlobalMembers.floor(ngsimd.GlobalMembers.log(dist[ord_dir[j]] / ngsimd.GlobalMembers.sqrt(2.0)) / ngsimd.GlobalMembers.log(fac1)),0.0);
		  }

		  // cout << " ref " << refi[0] << "\t" << refi[1] << "\t" << refi[2] << endl;
		  // cout << " order " << act_ref +1 - refi[0] << "\t" << act_ref +1 - refi[1] << "\t" << act_ref +1 - refi[2] << endl;

		  if (setorders)
		  {
			mesh[i].SetOrder(act_ref + 1 - refi[0],act_ref + 1 - refi[1],act_ref + 1 - refi[2]);
		  }
	  }
	  for (SurfaceElementIndex i = 0;i < mesh.GetNSE(); i++)
	  {
		  // Element2d el = mesh[i] ;
		  HPRefElement hpel = hpelements[mesh[i].hp_elnr];
		  ELEMENT_EDGE edges = MeshTopology.GetEdges1(mesh[i].GetType());
		  double[] dist = {0, 0, 0};
		  int[] ord_dir = {0, 0, 0};
		  int[] edge_dir = {0, 0, 0, 0};
		  int ned = 3;

		  if (mesh[i].GetType() == QUAD)
		  {
		  /*	cout << " QUAD " ;
		  for(int k=0;k<4;k++) cout << el[k] << "\t" ;
		  cout << endl; 	*/

		  edge_dir[2] = edge_dir[3] = 1;
		  ord_dir[1] = 1;
		  ned = 4;
		  }
		  /*  else
		    {
		  cout << " TRIG " ;
		  for(int k=0;k<3;k++) cout << el[k] << "\t" ;
		  cout << endl;
		  } */

		  for (int j = 0;j < ned;j++)
		  {
		  Vec < 3> v(hpel.param[edges[j][0] - 1][0] - hpel.param[edges[j][1] - 1][0], hpel.param[edges[j][0] - 1][1] - hpel.param[edges[j][1] - 1][1], hpel.param[edges[j][0] - 1][2] - hpel.param[edges[j][1] - 1][2]);
		  dist[edge_dir[j]] = Math.Max(v.Length(),dist[edge_dir[j]]);
		  }

		  int[] refi = new int[3];
		  for (int j = 0;j < 3;j++)
		  {
			refi[j] = (int)Math.Max((double)ngsimd.GlobalMembers.floor(ngsimd.GlobalMembers.log(dist[ord_dir[j]] / ngsimd.GlobalMembers.sqrt(2.0)) / ngsimd.GlobalMembers.log(fac1)),0.0);
		  }

		  if (setorders)
		  {
			mesh[i].SetOrder(act_ref + 1 - refi[0],act_ref + 1 - refi[1],act_ref + 1 - refi[2]);
		  }

			// cout << " ref " << refi[0] << "\t" << refi[1] << endl;
			// cout << " order " << act_ref +1 - refi[0] << "\t" << act_ref +1 - refi[1] << endl;
	  }
	  }
	}


//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void CalcTriangleBadness(double x2, double x3, double y3, double metricweight, double h, ref double badness, ref double g1x, ref double g1y);




//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//double CalcTriangleBadness(Point<3> p1, Point<3> p2, Point<3> p3, double metricweight, double h);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//double CalcTriangleBadness(Point<3> p1, Point<3> p2, Point<3> p3, Vec<3> n, double metricweight, double h);





	#if ! SMALLLIB
	//#ifndef NOTCL
	//#include <visual.hpp>
	//#endif
	#endif

	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_EXPORT __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_IMPORT __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_EXPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_IMPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max) { if ((value)<(min) || (value)>=(max)) throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max)); }
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_START(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_END(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_TRACER(n)


//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//double CalcTotalBad(Mesh::T_POINTS points, Array<Element, 0, uint> elements, MeshingParameters mp);


	public static double CalcBad(Mesh.T_POINTS points, Element elem, double h, MeshingParameters mp)
	{
	  if (elem.GetType() == TET)
	  {
		return CalcTetBadness(points[elem[0]], points[elem[1]], points[elem[2]], points[elem[3]], h, mp);
	  }
	  return 0;
	}



//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//int WrongOrientation(Mesh::T_POINTS points, Element el);

// extern double teterrpow; 
	// class CSGeometry;

	/// Build tet-mesh
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	static Timer MeshVolume_t("MeshVolume");

	public static MESHING3_RESULT MeshVolume(MeshingParameters mp, Mesh mesh3d)
	{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//  static Timer t("MeshVolume");
	  RegionTimer reg = new RegionTimer(MeshVolume_t);

	   int oldne;
	   int meshed;

	   Array<INDEX_2> connectednodes = new Array<INDEX_2>();

	   if (!mesh3d.HasLocalHFunction())
	   {
		   mesh3d.CalcLocalH(mp.grading);
	   }

	   mesh3d.Compress();

	   //  mesh3d.PrintMemInfo (cout);

	   if (mp.checkoverlappingboundary)
	   {
		  if (mesh3d.CheckOverlappingBoundary())
		  {
			 throw new Exception("Stop meshing since boundary mesh is overlapping");
		  }
	   }

	   int nonconsist = 0;
	   for (int k = 1; k <= mesh3d.GetNDomains(); k++)
	   {
		 if (mp.only3D_domain_nr && mp.only3D_domain_nr != k)
		 {
	   continue;
		 }
		  PrintMessage(3, "Check subdomain ", k, " / ", mesh3d.GetNDomains());

		  mesh3d.FindOpenElements(k);

		  /*
		  bool res = mesh3d.CheckOverlappingBoundary();
		  if (res)
		  {
		  PrintError ("Surface is overlapping !!");
		  nonconsist = 1;
		  }
		  */

		  bool res = (mesh3d.CheckConsistentBoundary() != 0);
		  if (res)
		  {
			 PrintError("Surface mesh not consistent");
			 nonconsist = 1;
		  }
	   }

	   if (nonconsist != 0)
	   {
		  PrintError("Stop meshing since surface mesh not consistent");
		  throw new Exception("Stop meshing since surface mesh not consistent");
	   }

	   double globmaxh = mp.maxh;

	   for (int k = 1; k <= mesh3d.GetNDomains(); k++)
	   {
	   if (mp.only3D_domain_nr && mp.only3D_domain_nr != k)
	   {
		 continue;
	   }
	   if (multithread.terminate)
	   {
			 break;
	   }

	   PrintMessage(2, "");
	   PrintMessage(1, "Meshing subdomain ", k, " of ", mesh3d.GetNDomains());
	   (*testout) << "Meshing subdomain " << k << "\n";

	   mp.maxh = netgen.GlobalMembers.min2(globmaxh, mesh3d.MaxHDomain(k));

	   mesh3d.CalcSurfacesOfNode();
	   mesh3d.FindOpenElements(k);

	   if (!mesh3d.GetNOpenElements())
	   {
			 continue;
	   }



	   Box < 3> domain_bbox(Box < 3>.EMPTY_BOX);

	   for (SurfaceElementIndex sei = 0; sei < mesh3d.GetNSE(); sei++)
	   {
		   Element2d el = mesh3d[sei];
		   if (el.IsDeleted())
		   {
			   continue;
		   }

		   if (mesh3d.GetFaceDescriptor(el.GetIndex()).DomainIn() == k || mesh3d.GetFaceDescriptor(el.GetIndex()).DomainOut() == k)
		   {

			 for (int j = 0; j < el.GetNP(); j++)
			 {
		   domain_bbox.Add(mesh3d[el[j]]);
			 }
		   }
	   }
	   domain_bbox.Increase(0.01 * domain_bbox.Diam());


		   for (int qstep = 0; qstep <= 3; qstep++)
		   {
			 // for (int qstep = 0; qstep <= 0; qstep++)  // for hex-filling
			  if (qstep == 0 && !mp.try_hexes)
			  {
				  continue;
			  }

		  // cout << "openquads = " << mesh3d.HasOpenQuads() << endl;
		  if (mesh3d.HasOpenQuads())
		  {
		  string rulefile = ngdir;

		  string[] rulep = null;
		  switch (qstep)
		  {
			case 0:
					  rulefile = "/Users/joachim/gitlab/netgen/rules/hexa.rls";
					  rulep = hexrules;
			  break;
			case 1:
			  rulefile += "/rules/prisms2.rls";
			  rulep = prismrules2;
			  break;
			case 2: // connect pyramid to triangle
			  rulefile += "/rules/pyramids2.rls";
			  rulep = pyramidrules2;
			  break;
			case 3: // connect to vis-a-vis point
			  rulefile += "/rules/pyramids.rls";
			  rulep = pyramidrules;
			  break;
		  }

				  // Meshing3 meshing(rulefile);
				  Meshing3 meshing = new Meshing3(rulep);

		  MeshingParameters mpquad = new MeshingParameters(mp);

		  mpquad.giveuptol = 15;
		  mpquad.baseelnp = 4;
		  mpquad.starshapeclass = 1000;
		  mpquad.check_impossible = qstep == 1; // for prisms only (air domain in trafo)


		  for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
		  {
			meshing.AddPoint(mesh3d[pi], pi);
		  }

				  /*
		  mesh3d.GetIdentifications().GetPairs (0, connectednodes);
		  for (int i = 1; i <= connectednodes.Size(); i++)
			meshing.AddConnectedPair (connectednodes.Get(i));
				  */
				  for (int nr = 1; nr <= mesh3d.GetIdentifications().GetMaxNr(); nr++)
				  {
					if (mesh3d.GetIdentifications().GetType(nr) != Identifications.PERIODIC)
					{
						mesh3d.GetIdentifications().GetPairs(nr, connectednodes);
						foreach (var pair in connectednodes)
						{
						  meshing.AddConnectedPair(pair);
						}
					}
				  }

		  for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
		  {
			  Element2d hel = mesh3d.OpenElement(i);
			  meshing.AddBoundaryElement(hel);
		  }

		  oldne = mesh3d.GetNE();

		  meshing.GenerateMesh(mesh3d, mpquad);

		  for (int i = oldne + 1; i <= mesh3d.GetNE(); i++)
		  {
			mesh3d.VolumeElement(i).SetIndex(k);
		  }

		  (*testout) << "mesh has " << mesh3d.GetNE() << " prism/pyramid elements" << "\n";

		  mesh3d.FindOpenElements(k);
		  }
		   }


		  if (mesh3d.HasOpenQuads())
		  {
			 PrintSysError("mesh has still open quads");
			 throw new Exception("Stop meshing since too many attempts");
			 // return MESHING3_GIVEUP;
		  }


		  if (mp.delaunay && mesh3d.GetNOpenElements())
		  {
			 Meshing3 meshing = new Meshing3((string)null);

			 mesh3d.FindOpenElements(k);

			 /*
			 for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
			    meshing.AddPoint (mesh3d[pi], pi);
			 */
			 foreach (PointIndex pi in mesh3d.Points().Range())
			 {
				meshing.AddPoint(mesh3d[pi], pi);
			 }

			 for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
			 {
				meshing.AddBoundaryElement(mesh3d.OpenElement(i));
			 }

			 oldne = mesh3d.GetNE();

			 meshing.Delaunay(mesh3d, k, mp);

			 for (int i = oldne + 1; i <= mesh3d.GetNE(); i++)
			 {
				mesh3d.VolumeElement(i).SetIndex(k);
			 }

			 PrintMessage(3, mesh3d.GetNP(), " points, ", mesh3d.GetNE(), " elements");
		  }


		  int cntsteps = 0;
		  if (mesh3d.GetNOpenElements())
		  {
			do
			{
				if (multithread.terminate)
				{
				  break;
				}

				mesh3d.FindOpenElements(k);
				PrintMessage(5, mesh3d.GetNOpenElements(), " open faces");
				cntsteps++;

				if (cntsteps > mp.maxoutersteps)
				{
				   throw new Exception("Stop meshing since too many attempts");
				}

				string rulefile = ngdir + "/tetra.rls";
				PrintMessage(1, "start tetmeshing");

				//	  Meshing3 meshing(rulefile);
				Meshing3 meshing = new Meshing3(tetrules);

				Array<int, PointIndex.BASE> glob2loc = new Array<int, PointIndex.BASE>(mesh3d.GetNP());
				glob2loc = -1;

				for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
				{
				  if (domain_bbox.IsIn(mesh3d[pi]))
				  {
					  glob2loc[pi] = meshing.AddPoint(mesh3d[pi], pi);
				  }
				}

				for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
				{
				   Element2d hel = mesh3d.OpenElement(i);
				   for (int j = 0; j < hel.GetNP(); j++)
				   {
					  hel[j] = glob2loc[hel[j]];
				   }
				   meshing.AddBoundaryElement(hel);
				   // meshing.AddBoundaryElement (mesh3d.OpenElement(i));
				}

				oldne = mesh3d.GetNE();

				mp.giveuptol = 15 + 10 * cntsteps;
				mp.sloppy = 5;
				meshing.GenerateMesh(mesh3d, mp);

				for (ElementIndex ei = oldne; ei < mesh3d.GetNE(); ei++)
				{
				   mesh3d[ei].SetIndex(k);
				}


				mesh3d.CalcSurfacesOfNode();
				mesh3d.FindOpenElements(k);

				// teterrpow = 2;
				if (mesh3d.GetNOpenElements() != 0)
				{
				   meshed = 0;
				   PrintMessage(5, mesh3d.GetNOpenElements(), " open faces found");

				   MeshOptimize3d optmesh = new MeshOptimize3d(mp);

				   const string optstr = "mcmstmcmstmcmstmcm";
				   for (uint j = 1; j <= optstr.Length; j++)
				   {
					  mesh3d.CalcSurfacesOfNode();
					  mesh3d.FreeOpenElementsEnvironment(2);
					  mesh3d.CalcSurfacesOfNode();

					  switch (optstr[j - 1])
					  {
					  case 'c':
						  optmesh.CombineImprove(mesh3d, OPT_REST);
						  break;
					  case 'd':
						  optmesh.SplitImprove(mesh3d, OPT_REST);
						  break;
					  case 's':
						  optmesh.SwapImprove(mesh3d, OPT_REST);
						  break;
					  case 't':
						  optmesh.SwapImprove2(mesh3d, OPT_REST);
						  break;
					  case 'm':
						  mesh3d.ImproveMesh(mp, OPT_REST);
						  break;
					  }

				   }

				   mesh3d.FindOpenElements(k);
				   PrintMessage(3, "Call remove problem");
				   RemoveProblem(mesh3d, k);
				   mesh3d.FindOpenElements(k);
				}
				else
				{
				   meshed = 1;
				   PrintMessage(1, "Success !");
				}
			} while (meshed == 0);
		  }

		  PrintMessage(1, mesh3d.GetNP(), " points, ", mesh3d.GetNE(), " elements");
	   }

	   mp.maxh = globmaxh;

	   MeshQuality3d(mesh3d);

	   return MESHING3_OK;
	}

/*


MESHING3_RESULT MeshVolumeOld (MeshingParameters & mp, Mesh& mesh3d)
{
int i, k, oldne;


int meshed;
int cntsteps; 


PlotStatistics3d * pstat;
if (globflags.GetNumFlag("silentflag", 1) <= 2)
pstat = new XPlotStatistics3d;
else
pstat = new TerminalPlotStatistics3d;

cntsteps = 0;
do
{
cntsteps++;
if (cntsteps > mp.maxoutersteps) 
{
return MESHING3_OUTERSTEPSEXCEEDED;
}


int noldp = mesh3d.GetNP();
    
    
if ( (cntsteps == 1) && globflags.GetDefineFlag ("delaunay"))
{
cntsteps ++;

mesh3d.CalcSurfacesOfNode();


for (k = 1; k <= mesh3d.GetNDomains(); k++)
{
Meshing3 meshing(NULL, pstat);

mesh3d.FindOpenElements(k);

for (i = 1; i <= noldp; i++)
meshing.AddPoint (mesh3d.Point(i), i);

for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
{
if (mesh3d.OpenElement(i).GetIndex() == k)
meshing.AddBoundaryElement (mesh3d.OpenElement(i));
}

oldne = mesh3d.GetNE();
if (globflags.GetDefineFlag ("blockfill"))
{
if (!globflags.GetDefineFlag ("localh"))
meshing.BlockFill 
(mesh3d, mp.h * globflags.GetNumFlag ("relblockfillh", 1));
else
meshing.BlockFillLocalH (mesh3d);
}

MeshingParameters mpd;
meshing.Delaunay (mesh3d, mpd);

for (i = oldne + 1; i <= mesh3d.GetNE(); i++)
mesh3d.VolumeElement(i).SetIndex (k);
}
}

noldp = mesh3d.GetNP();

mesh3d.CalcSurfacesOfNode();
mesh3d.FindOpenElements();
for (k = 1; k <= mesh3d.GetNDomains(); k++)
{
Meshing3 meshing(globflags.GetStringFlag ("rules3d", NULL), pstat);
    
Point3d pmin, pmax;
mesh3d.GetBox (pmin, pmax, k);

rot.SetCenter (Center (pmin, pmax));

for (i = 1; i <= noldp; i++)
meshing.AddPoint (mesh3d.Point(i), i);

for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
{
if (mesh3d.OpenElement(i).GetIndex() == k)
meshing.AddBoundaryElement (mesh3d.OpenElement(i));
}

oldne = mesh3d.GetNE();


if ( (cntsteps == 1) && globflags.GetDefineFlag ("blockfill"))
{
if (!globflags.GetDefineFlag ("localh"))
{
meshing.BlockFill 
(mesh3d, 
mp.h * globflags.GetNumFlag ("relblockfillh", 1));
}
else
{
meshing.BlockFillLocalH (mesh3d);
}
}


mp.giveuptol = int(globflags.GetNumFlag ("giveuptol", 15));

meshing.GenerateMesh (mesh3d, mp);

for (i = oldne + 1; i <= mesh3d.GetNE(); i++)
mesh3d.VolumeElement(i).SetIndex (k);
}



mesh3d.CalcSurfacesOfNode();
mesh3d.FindOpenElements();
    
teterrpow = 2;
if (mesh3d.GetNOpenElements() != 0)
{
meshed = 0;
(*mycout) << "Open elements found, old" << endl;
const char * optstr = "mcmcmcmcm";
int j;
for (j = 1; j <= strlen(optstr); j++)
switch (optstr[j-1])
{
case 'c': mesh3d.CombineImprove(); break;
case 'd': mesh3d.SplitImprove(); break;
case 's': mesh3d.SwapImprove(); break;
case 'm': mesh3d.ImproveMesh(2); break;
}	  

(*mycout) << "Call remove" << endl;
RemoveProblem (mesh3d);
(*mycout) << "Problem removed" << endl;
}
else
meshed = 1;
}
while (!meshed);

MeshQuality3d (mesh3d);

return MESHING3_OK;
}  

*/




/*
MESHING3_RESULT MeshMixedVolume(MeshingParameters & mp, Mesh& mesh3d)
{
  int i, j;
  MESHING3_RESULT res;
  Point3d pmin, pmax;

  mp.giveuptol = 10;
  mp.baseelnp = 4;
  mp.starshapeclass = 100;

  //  TerminalPlotStatistics3d pstat;

  Meshing3 meshing1("pyramids.rls");
  for (i = 1; i <= mesh3d.GetNP(); i++)
    meshing1.AddPoint (mesh3d.Point(i), i);

  mesh3d.FindOpenElements();
  for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
    if (mesh3d.OpenElement(i).GetIndex() == 1)
  meshing1.AddBoundaryElement (mesh3d.OpenElement(i));

  res = meshing1.GenerateMesh (mesh3d, mp);

  mesh3d.GetBox (pmin, pmax);
  PrintMessage (1, "Mesh pyramids, res = ", res);
  if (res)
    exit (1);


  for (i = 1; i <= mesh3d.GetNE(); i++)
    mesh3d.VolumeElement(i).SetIndex (1);

  // do delaunay

  mp.baseelnp = 0;
  mp.starshapeclass = 5;

  Meshing3 meshing2(NULL);
  for (i = 1; i <= mesh3d.GetNP(); i++)
    meshing2.AddPoint (mesh3d.Point(i), i);
  
  mesh3d.FindOpenElements();
  for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
    if (mesh3d.OpenElement(i).GetIndex() == 1)
  meshing2.AddBoundaryElement (mesh3d.OpenElement(i));

  MeshingParameters mpd;
  meshing2.Delaunay (mesh3d, mpd);

  for (i = 1; i <= mesh3d.GetNE(); i++)
    mesh3d.VolumeElement(i).SetIndex (1);


  mp.baseelnp = 0;
  mp.giveuptol = 10;

  for (int trials = 1; trials <= 50; trials++)
    {
  if (multithread.terminate)
	return MESHING3_TERMINATE;

  Meshing3 meshing3("tetra.rls");
  for (i = 1; i <= mesh3d.GetNP(); i++)
	meshing3.AddPoint (mesh3d.Point(i), i);
    
  mesh3d.FindOpenElements();
  for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
	if (mesh3d.OpenElement(i).GetIndex() == 1)
	  meshing3.AddBoundaryElement (mesh3d.OpenElement(i));
    
  if (trials > 1)
	CheckSurfaceMesh2 (mesh3d);
  res = meshing3.GenerateMesh (mesh3d, mp);
    
  for (i = 1; i <= mesh3d.GetNE(); i++)
	mesh3d.VolumeElement(i).SetIndex (1);

  if (res == 0) break;



  for (i = 1; i <= mesh3d.GetNE(); i++)
	{
	  const Element & el = mesh3d.VolumeElement(i);
	  if (el.GetNP() != 4)
		{
	  for (j = 1; j <= el.GetNP(); j++)
		mesh3d.AddLockedPoint (el.PNum(j));
		}
	}

  mesh3d.CalcSurfacesOfNode();
  mesh3d.FindOpenElements();

  MeshOptimize3d optmesh;

  teterrpow = 2;
  const char * optstr = "mcmcmcmcm";
  for (j = 1; j <= strlen(optstr); j++)
	switch (optstr[j-1])
	  {
	  case 'c': optmesh.CombineImprove(mesh3d, OPT_REST); break;
	  case 'd': optmesh.SplitImprove(mesh3d); break;
	  case 's': optmesh.SwapImprove(mesh3d); break;
	  case 'm': mesh3d.ImproveMesh(); break;
	  }	  

  RemoveProblem (mesh3d);
    }


  PrintMessage (1, "Meshing tets, res = ", res);
  if (res)
    {
  mesh3d.FindOpenElements();
  PrintSysError (1, "Open elements: ", mesh3d.GetNOpenElements());
  exit (1);
    }



  for (i = 1; i <= mesh3d.GetNE(); i++)
    {
  const Element & el = mesh3d.VolumeElement(i);
  if (el.GetNP() != 4)
	{
	  for (j = 1; j <= el.GetNP(); j++)
		mesh3d.AddLockedPoint (el.PNum(j));
	}
    }

  mesh3d.CalcSurfacesOfNode();
  mesh3d.FindOpenElements();

  MeshOptimize3d optmesh;

  teterrpow = 2;
  const char * optstr = "mcmcmcmcm";
  for (j = 1; j <= strlen(optstr); j++)
    switch (optstr[j-1])
  {
  case 'c': optmesh.CombineImprove(mesh3d, OPT_REST); break;
  case 'd': optmesh.SplitImprove(mesh3d); break;
  case 's': optmesh.SwapImprove(mesh3d); break;
  case 'm': mesh3d.ImproveMesh(); break;
  }	  


  return MESHING3_OK;
}
*/







	/// Build mixed-element mesh
	// MESHING3_RESULT MeshMixedVolume (MeshingParameters & mp, Mesh& mesh3d);

	/// Optimize tet-mesh
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	static Timer OptimizeVolume_t("OptimizeVolume");

	public static MESHING3_RESULT OptimizeVolume(MeshingParameters mp, Mesh mesh3d)
	{
	  //				  const CSGeometry * geometry)
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//  static Timer t("OptimizeVolume");
	  RegionTimer reg = new RegionTimer(OptimizeVolume_t);

	  int i;

	  PrintMessage(1, "Volume Optimization");

	  /*
	    if (!mesh3d.PureTetMesh())
	    return MESHING3_OK;
	  */

	  // (*mycout) << "optstring = " << mp.optimize3d << endl;
	  /*
	    const char * optstr = globflags.GetStringFlag ("optimize3d", "cmh");
	    int optsteps = int (globflags.GetNumFlag ("optsteps3d", 2));
	  */

	  mesh3d.CalcSurfacesOfNode();
	  for (i = 1; i <= mp.optsteps3d; i++)
	  {
	  if (multithread.terminate)
	  {
		break;
	  }

	  MeshOptimize3d optmesh = new MeshOptimize3d(mp);

	  // teterrpow = mp.opterrpow;
	  // for (size_t j = 1; j <= strlen(mp.optimize3d); j++)
		  for (uint j = 1; j <= mp.optimize3d.length(); j++)
		  {
		  if (multithread.terminate)
		  {
			break;
		  }

		  switch (mp.optimize3d[j - 1])
		  {
			case 'c':
				optmesh.CombineImprove(mesh3d, OPT_REST);
				break;
			case 'd':
				optmesh.SplitImprove(mesh3d);
				break;
			case 's':
				optmesh.SwapImprove(mesh3d);
				break;
				  // case 'u': optmesh.SwapImproveSurface(mesh3d); break;
			case 't':
				optmesh.SwapImprove2(mesh3d);
				break;
	#if SOLIDGEOM
//C++ TO C# CONVERTER TODO TASK: C# does not allow fall-through from a non-empty 'case':
			case 'm':
				mesh3d.ImproveMesh(*geometry);
				break;
			case 'M':
				mesh3d.ImproveMesh(*geometry);
				break;
	#else
//C++ TO C# CONVERTER TODO TASK: C# does not allow fall-through from a non-empty 'case':
			case 'm':
				mesh3d.ImproveMesh(mp);
				break;
			case 'M':
				mesh3d.ImproveMesh(mp);
				break;
	#endif
//C++ TO C# CONVERTER TODO TASK: C# does not allow fall-through from a non-empty 'case':
			case 'j':
				mesh3d.ImproveMeshJacobian(mp);
				break;
		  }
		  }
	  mesh3d.mglevels = 1;
	  MeshQuality3d(mesh3d);
	  }

	  return MESHING3_OK;
	}

	//			       const CSGeometry * geometry = NULL);

	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	static Timer RemoveIllegalElements_t("RemoveIllegalElements");

	public static void RemoveIllegalElements(Mesh mesh3d)
	{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//  static Timer t("RemoveIllegalElements");
	  RegionTimer reg = new RegionTimer(RemoveIllegalElements_t);

	  int it = 10;
	  int nillegal;
	  int oldn;

	  PrintMessage(1, "Remove Illegal Elements");
	  // return, if non-pure tet-mesh
	  /*
	    if (!mesh3d.PureTetMesh())
	    return;
	  */
	  mesh3d.CalcSurfacesOfNode();

	  nillegal = mesh3d.MarkIllegalElements();

	  MeshingParameters dummymp = new MeshingParameters();
	  MeshOptimize3d optmesh = new MeshOptimize3d(dummymp);
	  while (nillegal != 0 && (it--) > 0)
	  {
	  if (multithread.terminate)
	  {
		break;
	  }

	  PrintMessage(5, nillegal, " illegal tets");
		  optmesh.SplitImprove(mesh3d, OPT_LEGAL);

	  mesh3d.MarkIllegalElements(); // test
	  optmesh.SwapImprove(mesh3d, OPT_LEGAL);
	  mesh3d.MarkIllegalElements(); // test
	  optmesh.SwapImprove2(mesh3d, OPT_LEGAL);

	  oldn = nillegal;
	  nillegal = mesh3d.MarkIllegalElements();

	  if (oldn != nillegal)
	  {
		it = 10;
	  }
	  }
	  PrintMessage(5, nillegal, " illegal tets");
	}

	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_EXPORT __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_IMPORT __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_EXPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_IMPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max) { if ((value)<(min) || (value)>=(max)) throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max)); }
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_START(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_END(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_TRACER(n)


	///
	public static void MeshQuality2d(Mesh mesh)
	{
	  int ncl = 20;
	  int cl;
	  Array<int> incl = new Array<int>(ncl);
	  int i;
	  SurfaceElementIndex sei = new SurfaceElementIndex();
	  double qual;

	  incl = 0;

	  for (sei = 0; sei < mesh.GetNSE(); sei++)
	  {
	  qual = netgen.GlobalMembers.TriangleQualityInst(mesh[mesh[sei][0]], mesh[mesh[sei][1]], mesh[mesh[sei][2]]);

	  cl = (int)((ncl - 1e-3) * qual) + 1;
	  incl.Elem(cl)++;
	  }

	  (*testout) << "\n" << "\n";

	  (*testout) << "Points:           " << mesh.GetNP() << "\n";
	  (*testout) << "Surface Elements: " << mesh.GetNSE() << "\n";

	  (*testout) << "\n";
	  (*testout) << "Elements in qualityclasses:" << "\n";
	  // (*testout).precision(2);
	  (*testout) << setprecision(2);
	  for (i = 1; i <= ncl; i++)
	  {
	  (*testout) << setw(4) << (double)(i - 1) / ncl << " - " << setw(4) << (double)i / ncl << ": " << incl.Get(i) << "\n";
	  }
	}

	///
	public static void MeshQuality3d(Mesh mesh, Array<int> inclass = null)
	{
	  int ncl = 20;
	  int cl;
	  Array<int> incl = new Array<int>(ncl);
	  int i;
	  double qual;
	  double sum = 0;
	  int nontet = 0;

	  for (i = 1; i <= incl.Size(); i++)
	  {
		incl.Elem(i) = 0;
	  }

	  for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	  {
	  if (mesh[ei].GetType() != TET)
	  {
		  nontet++;
		  continue;
	  }

	  qual = netgen.GlobalMembers.TetElementQuality(new mesh.Point(mesh[ei][0]), new mesh.Point(mesh[ei][1]), new mesh.Point(mesh[ei][2]), new mesh.Point(mesh[ei][3]));

	  if (qual > 1)
	  {
		  qual = 1;
	  }
	  cl = (int)(ncl * qual) + 1;

	  if (cl < 1)
	  {
		  cl = 1;
	  }
	  if (cl > ncl)
	  {
		  cl = ncl;
	  }

	  incl.Elem(cl)++;
	  if (inclass != null)
	  {
		  inclass[ei] = cl;
	  }
	  sum += 1 / qual;
	  }

	  (*testout) << "\n" << "\n";
	  (*testout) << "Points:           " << mesh.GetNP() << "\n";
	  (*testout) << "Volume Elements:  " << mesh.GetNE() << "\n";
	  if (nontet != 0)
	  {
		(*testout) << nontet << " non tetrahedral elements" << "\n";
	  }
	  (*testout) << "\n";

	  (*testout) << "Volume elements in qualityclasses:" << "\n";
	  (*testout) << setprecision(2);
	  for (i = 1; i <= ncl; i++)
	  {
	  (*testout) << setw(4) << (double)(i - 1) / ncl << " - " << setw(4) << (double)i / ncl << ": " << incl.Get(i) << "\n";
	  }
	  (*testout) << "total error: " << sum << "\n";
	}

	///
	public static void SaveEdges(Mesh mesh, string geomfile, double h, ref string filename)
	{
	  ofstream of = new ofstream(filename);
	  int i;
	  Segment seg;

	  of << "edges" << "\n";
	  of << geomfile << "\n";
	  of << h << "\n";

	  of << mesh.GetNP() << "\n";
	  for (i = 1; i <= mesh.GetNP(); i++)
	  {
		of << new mesh.Point(i)(0) << " " << new mesh.Point(i)(1) << " " << new mesh.Point(i)(2) << "\n";
	  }

	  of << 2 * mesh.GetNSeg() << "\n";
	  for (i = 1; i <= mesh.GetNSeg(); i++)
	  {
	  seg = mesh.LineSegment(i);

	  of << seg[1] << " " << seg[0] << " " << seg.si << "\n";
	  }

	}

	///
	public static void SaveSurfaceMesh(Mesh mesh, double h, ref string filename)

	{
	  int i;

	  ofstream outfile = new ofstream(filename);

	  outfile << "surfacemesh" << "\n";
	  outfile << h << "\n";

	  outfile << mesh.GetNP() << "\n";
	  for (i = 1; i <= mesh.GetNP(); i++)
	  {
		outfile << new mesh.Point(i)(0) << " " << new mesh.Point(i)(1) << " " << new mesh.Point(i)(2) << "\n";
	  }



	  outfile << mesh.GetNSE() << "\n";
	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	  Element2d el = mesh.SurfaceElement(i);

	  if (mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
	  {
		outfile << mesh.SurfaceElement(i).PNum(1) << " " << mesh.SurfaceElement(i).PNum(2) << " " << mesh.SurfaceElement(i).PNum(3) << "\n";
	  }
	  if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0)
	  {
		outfile << mesh.SurfaceElement(i).PNum(1) << " " << mesh.SurfaceElement(i).PNum(3) << " " << mesh.SurfaceElement(i).PNum(2) << "\n";
	  }
	  }
	}

	/*
	///
	extern void Save2DMesh (
	         const Mesh & mesh2d,
		 const Array<class SplineSegment*> * splines,
	         ostream & outfile);
	*/

//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
	//class Surface;
	///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void SaveVolumeMesh(Array<Point3d> points, Array<Element> elements, Array<Element> volelements, Array<Surface> surfaces, ref string filename);

	///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void SaveVolumeMesh(Mesh mesh, NetgenGeometry geometry, ref string filename);

	///
	public static int CheckCode()
	{
	  return 1;

	  /*
	    char st[100];
	    ifstream ist("pw");
	
	    if (!ist.good()) return 0;
	    ist >> st;
	    if (strcmp (st, "JKULinz") == 0) return 1;
	    return 0;
	  */
	}

// static double teterrpow = 2;



	///
	public static double CalcTetBadness(Point3d p1, Point3d p2, Point3d p3, Point3d p4, double h, MeshingParameters mp)
	{
	  double vol;
	  double l;
	  double ll;
	  double lll;
	  double ll1;
	  double ll2;
	  double ll3;
	  double ll4;
	  double ll5;
	  double ll6;
	  double err;

	  Vec3d v1 = new Vec3d(p1, p2);
	  Vec3d v2 = new Vec3d(p1, p3);
	  Vec3d v3 = new Vec3d(p1, p4);

	  vol = netgen.GlobalMembers.Determinant(v1, v2, v3) * (-0.166666666666666);

	  ll1 = v1.Length2();
	  ll2 = v2.Length2();
	  ll3 = v3.Length2();
	  ll4 = netgen.GlobalMembers.Dist2(p2, p3);
	  ll5 = netgen.GlobalMembers.Dist2(p2, p4);
	  ll6 = netgen.GlobalMembers.Dist2(p3, p4);

	  ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
	  l = ngsimd.GlobalMembers.sqrt(ll);
	  lll = l * ll;

	  if (vol <= 1e-24 * lll)
	  {
		return 1e24;
	  }

	  err = 0.0080187537 * lll / vol; // sqrt(216) / (6^4 * sqrt(2))

	  if (h > 0)
	  {
		err += ll / (h * h) + h * h * (1 / ll1 + 1 / ll2 + 1 / ll3 + 1 / ll4 + 1 / ll5 + 1 / ll6) - 12;
	  }

	  double teterrpow = mp.opterrpow;
	  if (teterrpow < 1)
	  {
		  teterrpow = 1;
	  }

	  if (teterrpow == 1)
	  {
		  return err;
	  }
	  if (teterrpow == 2)
	  {
		  return err * err;
	  }
	  return new SIMD<double,N>(ngsimd.GlobalMembers.pow(err, teterrpow));
	}

	///
	public static double CalcTetBadnessGrad(Point3d p1, Point3d p2, Point3d p3, Point3d p4, double h, int pi, ref Vec < 3> grad, MeshingParameters mp)
	{
	  double vol;
	  double l;
	  double ll;
	  double lll;
	  double err;

	  Point3d pp1;
	  Point3d pp2;
	  Point3d pp3;
	  Point3d pp4;

	  pp1 = p1;
	  pp2 = p2;
	  pp3 = p3;
	  pp4 = p4;

	  switch (pi)
	  {
		case 2:
		{
		swap(pp1, pp2);
		swap(pp3, pp4);
		break;
		}
		case 3:
		{
		swap(pp1, pp3);
		swap(pp2, pp4);
		break;
		}
		case 4:
		{
		swap(pp1, pp4);
		swap(pp3, pp2);
		break;
		}
	  }


	  Vec3d v1 = new Vec3d(pp1, pp2);
	  Vec3d v2 = new Vec3d(pp1, pp3);
	  Vec3d v3 = new Vec3d(pp1, pp4);

	  Vec3d v4 = new Vec3d(pp2, pp3);
	  Vec3d v5 = new Vec3d(pp2, pp4);
	  Vec3d v6 = new Vec3d(pp3, pp4);

	  vol = netgen.GlobalMembers.Determinant(v1, v2, v3) * (-0.166666666666666);

	  Vec3d gradvol = new Vec3d();
	  netgen.GlobalMembers.Cross(v5, v4, gradvol);
	  gradvol *= (-1.0 / 6.0);


	  double ll1 = v1.Length2();
	  double ll2 = v2.Length2();
	  double ll3 = v3.Length2();
	  double ll4 = v4.Length2();
	  double ll5 = v5.Length2();
	  double ll6 = v6.Length2();

	  ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
	  l = ngsimd.GlobalMembers.sqrt(ll);
	  lll = l * ll;

	  if (vol <= 1e-24 * lll)
	  {
	  grad = new Vec3d(0, 0, 0);
	  return 1e24;
	  }



	  Vec3d gradll1 = new Vec3d(pp2, pp1);
	  Vec3d gradll2 = new Vec3d(pp3, pp1);
	  Vec3d gradll3 = new Vec3d(pp4, pp1);
	  gradll1 *= 2;
	  gradll2 *= 2;
	  gradll3 *= 2;

	  Vec3d gradll = new Vec3d(gradll1);
	  gradll += gradll2;
	  gradll += gradll3;

	  /*
	  Vec3d gradll;
	  gradll = v1+v2+v3;
	  gradll *= -2;
	  */

	  err = 0.0080187537 * lll / vol;


	  gradll *= (0.0080187537 * 1.5 * l / vol);
	  Vec3d graderr = new Vec3d(gradll);
	  gradvol *= (-0.0080187537 * lll / (vol * vol));
	  graderr += gradvol;

	  if (h > 0)
	  {
	  /*
	  Vec3d gradll1 (*pp2, *pp1);
	  Vec3d gradll2 (*pp3, *pp1);
	  Vec3d gradll3 (*pp4, *pp1);
	  gradll1 *= 2;
	  gradll2 *= 2;
	  gradll3 *= 2;
	  */
	  err += ll / (h * h) + h * h * (1 / ll1 + 1 / ll2 + 1 / ll3 + 1 / ll4 + 1 / ll5 + 1 / ll6) - 12;

	  graderr += (1 / (h * h) - h * h / (ll1 * ll1)) * gradll1;
	  graderr += (1 / (h * h) - h * h / (ll2 * ll2)) * gradll2;
	  graderr += (1 / (h * h) - h * h / (ll3 * ll3)) * gradll3;
	  }

	  double errpow;

	  double teterrpow = mp.opterrpow;
	  if (teterrpow < 1)
	  {
		  teterrpow = 1;
	  }

	  if (teterrpow == 1)
	  {
		  errpow = err;
		  grad = graderr;
	  }
	  else if (teterrpow == 2)
	  {
		  errpow = err * err;
		  grad = (2 * err) * graderr;
	  }
	  else
	  {
		  errpow = ngsimd.GlobalMembers.pow(err, teterrpow);
		  grad = (teterrpow * errpow / err) * graderr;
	  }
	  return errpow;
	}

/*

double CalcTetBadness (const Point3d & p1, const Point3d & p2,
const Point3d & p3, const Point3d & p4, double h)
{
double vol, l;
double err;


Vec3d v1 (p1, p2);
Vec3d v2 (p1, p3);
Vec3d v3 (p1, p4);

vol = -Determinant (v1, v2, v3) / 6;

double l1 = v1.Length();
double l2 = v2.Length();
double l3 = v3.Length();
double l4 = Dist (p2, p3);
double l5 = Dist (p2, p4);
double l6 = Dist (p3, p4);

l = l1 + l2 + l3 + l4 + l5 + l6;

// just for timing
// l += 1e-40 * CalcTetBadnessNew (p1, p2, p3, p4, h);

if (vol <= 1e-24 * l * l * l)
{ 
return 1e24;
}

err = (l*l*l) / (1832.82 * vol);    // 6^4 * sqrt(2)

if (h > 0)
err += l / h + 
h * (1 / l1 + 1/l2 + 1/l3 + 1/l4 + 1/l5 + 1/l6) - 12;

return pow (err, teterrpow);
}



double CalcTetBadnessGrad (const Point3d & p1, const Point3d & p2,
const Point3d & p3, const Point3d & p4, double h,
int pi, Vec3d & grad)
{
double vol, l;
double err;

const Point3d *pp1, *pp2, *pp3, *pp4;

pp1 = &p1;
pp2 = &p2;
pp3 = &p3;
pp4 = &p4;

switch (pi)
{
case 2:
{
swap (pp1, pp2);
swap (pp3, pp4);
break;
}
case 3:
{
swap (pp1, pp3);
swap (pp2, pp4);
break;
}
case 4:
{
swap (pp1, pp4);
swap (pp3, pp2);
break;
}
}


Vec3d v1 (*pp1, *pp2);
Vec3d v2 (*pp1, *pp3);
Vec3d v3 (*pp1, *pp4);

Vec3d v4 (*pp2, *pp3);
Vec3d v5 (*pp2, *pp4);
Vec3d v6 (*pp3, *pp4);


//   Vec3d n;
//   Cross (v1, v2, n);
//   vol = - (n * v3) / 6;


vol = -Determinant (v1, v2, v3) / 6;  

Vec3d gradvol;
Cross (v5, v4, gradvol);
gradvol *= (-1.0/6.0);


double l1 = v1.Length();
double l2 = v2.Length();
double l3 = v3.Length();
double l4 = v4.Length();
double l5 = v5.Length();
double l6 = v6.Length();

l = l1 + l2 + l3 +l4 + l5 + l6;

Vec3d gradl1 (*pp2, *pp1);
Vec3d gradl2 (*pp3, *pp1);
Vec3d gradl3 (*pp4, *pp1);
gradl1 /= l1;
gradl2 /= l2;
gradl3 /= l3;

Vec3d gradl (gradl1);
gradl += gradl2;
gradl += gradl3;


if (vol <= 1e-24 * l * l * l)
{ 
grad = Vec3d (0, 0, 0);
return 1e24;
}


double c1 = 1.0 / 1832.82;      // 6^4 * sqrt(2)
err = c1 * (l*l*l) / vol; 


gradl *= (c1 * 3 * l * l / vol);
Vec3d graderr(gradl);
gradvol *= ( -c1 * l * l * l / (vol * vol) );
graderr+= gradvol;

if (h > 0)
{
err += l / h + 
h * ( 1 / l1 + 1 / l2 + 1 / l3 + 
1 / l4 + 1 / l5 + 1 / l6 ) - 12;

graderr += (1/h - h/(l1*l1)) * gradl1;
graderr += (1/h - h/(l2*l2)) * gradl2;
graderr += (1/h - h/(l3*l3)) * gradl3;
cout << "?";
}

double errpow = pow (err, teterrpow);
grad = (teterrpow * errpow / err) * graderr;

return errpow;
}

*/





/*
  double CalcVolume (const Array<Point3d> & points,
  const Element & el)
  {
  Vec3d v1 = points.Get(el.PNum(2)) - 
  points.Get(el.PNum(1));
  Vec3d v2 = points.Get(el.PNum(3)) - 
  points.Get(el.PNum(1));
  Vec3d v3 = points.Get(el.PNum(4)) - 
  points.Get(el.PNum(1)); 
       
  return -(Cross (v1, v2) * v3) / 6;	 
  }  
*/



	/** Calculates volume of an element.
	  The volume of the tetrahedron el is computed
	 */
	// extern double CalcVolume (const Array<Point3d> & points,
	//        const Element & el);  

	/** The total volume of all elements is computed.
	  This function calculates the volume of the mesh */
	public static double CalcVolume(Array<Point3d> points, Array<Element> elements)
	{
	  double vol;
	  Vec3d v1 = new Vec3d();
	  Vec3d v2 = new Vec3d();
	  Vec3d v3 = new Vec3d();

	  vol = 0;
	  for (int i = 0; i < elements.Size(); i++)
	  {
	  v1 = points.Get(elements[i][1]) - points.Get(elements[i][0]);
	  v2 = points.Get(elements[i][2]) - points.Get(elements[i][0]);
	  v3 = points.Get(elements[i][3]) - points.Get(elements[i][0]);
	  vol -= (netgen.GlobalMembers.Cross(v1, v2) * v3) / 6;
	  }
	  return vol;
	}

	///
	public static int CheckSurfaceMesh(Mesh mesh)
	{
	  PrintMessage(3, "Check Surface mesh");

	  int nf = mesh.GetNSE();
	  INDEX_2_HASHTABLE<int> edges = new INDEX_2_HASHTABLE<int>(nf + 2);
	  int i;
	  int j;
	  INDEX_2 i2 = new INDEX_2();
	  int cnt1 = 0;
	  int cnt2 = 0;

	  for (i = 1; i <= nf; i++)
	  {
		for (j = 1; j <= 3; j++)
		{
		i2.I1() = mesh.SurfaceElement(i).PNumMod(j);
		i2.I2() = mesh.SurfaceElement(i).PNumMod(j + 1);
		if (edges.Used(i2))
		{
			int hi;
			hi = edges.Get(i2);
			if (hi != 1)
			{
		  PrintSysError("CheckSurfaceMesh, hi = ", hi);
			}
			edges.Set(i2, 2);
			cnt2++;
		}
		else
		{
			netgen.GlobalMembers.Swap(ref i2.I1(), ref i2.I2());
			edges.Set(i2, 1);
			cnt1++;
		}
		}
	  }


	  if (cnt1 != cnt2)
	  {
	  PrintUserError("Surface mesh not consistent");
	  //      MyBeep(2);
	  //      (*mycout) << "cnt1 = " << cnt1 << " cnt2 = " << cnt2 << endl;
	  return 0;
	  }
	  return 1;
	}

	///
	public static int CheckSurfaceMesh2(Mesh mesh)
	{
	  int i;
	  int j;
	  int k;
	  const Point < 3> *tri1[3], *tri2[3];

	  for (i = 1; i <= mesh.GetNOpenElements(); i++)
	  {
	  PrintDot();
	  for (j = 1; j < i; j++)
	  {
		  for (k = 1; k <= 3; k++)
		  {
		  tri1[k - 1] = &new mesh.Point(mesh.OpenElement(i).PNum(k));
		  tri2[k - 1] = &new mesh.Point(mesh.OpenElement(j).PNum(k));
		  }
		  if (IntersectTriangleTriangle(tri1[0], tri2[0]))
		  {
		  PrintSysError("Surface elements are intersecting");
		  (*testout) << "Intersecting: " << "\n";
		  for (k = 0; k <= 2; k++)
		  {
			(*testout) << *tri1[k] << "   ";
		  }
		  (*testout) << "\n";
		  for (k = 0; k <= 2; k++)
		  {
			(*testout) << *tri2[k] << "   ";
		  }
		  (*testout) << "\n";
		  }

	  }
	  }
	  return 0;
	}

/* ******************** CheckMesh ******************************* */

/// Checks, whether mesh contains a valid 3d mesh

	///
	public static int CheckMesh3D(Mesh mesh)
	{
	  INDEX_3_HASHTABLE<int> faceused = new INDEX_3_HASHTABLE<int>(mesh.GetNE() / 3);
	  int i;
	  int j;
	  int k;
	  int l;
	  INDEX_3 i3 = new INDEX_3();
	  int ok = 1;
	  ElementIndex ei = new ElementIndex();

	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	  Element2d el = mesh.SurfaceElement(i);

	  if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0 || mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
	  {
		  for (j = 1; j <= 3; j++)
		  {
			i3.I(j) = el.PNum(j);
		  }

		  i3.Sort();
		  faceused.Set(i3, 1);
	  }
	  }

	  for (ei = 0; ei < mesh.GetNE(); ei++)
	  {
	  Element el = mesh[ei];

	  for (j = 1; j <= 4; j++)
	  {
		  l = 0;
		  for (k = 1; k <= 4; k++)
		  {
		  if (j != k)
		  {
			  l++;
			  i3.I(l) = el.PNum(k);
		  }
		  }

		  i3.Sort();
		  if (faceused.Used(i3))
		  {
			faceused.Set(i3, faceused.Get(i3) + 1);
		  }
		  else
		  {
			faceused.Set(i3, 1);
		  }
	  }
	  }


	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	  Element2d el = mesh.SurfaceElement(i);

	  for (j = 1; j <= 3; j++)
	  {
		i3.I(j) = el.PNum(j);
	  }

	  i3.Sort();
	  k = faceused.Get(i3);
	  if (k != 2)
	  {
		  ok = 0;
		  (*testout) << "face " << i << " with points " << i3.I1() << "-" << i3.I2() << "-" << i3.I3() << " has " << k << " elements" << "\n";
	  }
	  }

	  for (ei = 0; ei < mesh.GetNE(); ei++)
	  {
	  Element el = mesh[ei];

	  for (j = 1; j <= 4; j++)
	  {
		  l = 0;
		  for (k = 1; k <= 4; k++)
		  {
		  if (j != k)
		  {
			  l++;
			  i3.I(l) = el.PNum(k);
		  }
		  }

		  i3.Sort();
		  k = faceused.Get(i3);
		  if (k != 2)
		  {
		  ok = 0;
		  (*testout) << "element " << ei << " with face " << i3.I1() << "-" << i3.I2() << "-" << i3.I3() << " has " << k << " elements" << "\n";
		  }
	  }
	  }





	  /*
	    for (i = 1; i <= faceused.GetNBags(); i++)
	    for (j = 1; j <= faceused.GetBagSize(i); j++)
	    {
	    faceused.GetData(i, j, i3, k);
	    if (k != 2)
	    {
	    (*testout) << "Face: " << i3.I1() << "-"
	    << i3.I2() << "-" << i3.I3() << " has "
	    << k << " Faces " << endl;
	    cerr << "Face Error" << endl;
	    ok = 0;
	    }
	    }
	  */


	  if (ok == 0)
	  {
	  (*testout) << "surfelements: " << "\n";
	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
		  Element2d el = mesh.SurfaceElement(i);
		  (*testout) << setw(5) << i << ":" << setw(6) << el.GetIndex() << setw(6) << el.PNum(1) << setw(4) << el.PNum(2) << setw(4) << el.PNum(3) << "\n";
	  }
	  (*testout) << "volelements: " << "\n";
	  for (ei = 0; ei < mesh.GetNE(); ei++)
	  {
		  Element el = mesh[ei];
		  (*testout) << setw(5) << i << ":" << setw(6) << el.GetIndex() << setw(6) << el[0] << setw(4) << el[1] << setw(4) << el[2] << setw(4) << el[3] << "\n";
	  }
	  }


	  return ok;
	}

	///
	public static void RemoveProblem(Mesh mesh, int domainnr)
	{
	  int i;
	  int j;
	  int k;

	  mesh.FindOpenElements(domainnr);
	  int np = mesh.GetNP();

	  BitArrayChar<PointIndex.BASE> ppoints = new BitArrayChar<PointIndex.BASE>(np);

	  // int ndom = mesh.GetNDomains();

	  PrintMessage(3, "Elements before Remove: ", mesh.GetNE());
	  // for (k = 1; k <= ndom; k++)
	  k = domainnr;
	  {
	  ppoints.Clear();

	  for (i = 1; i <= mesh.GetNOpenElements(); i++)
	  {
		  Element2d sel = mesh.OpenElement(i);
		  if (sel.GetIndex() == k)
		  {
		  for (j = 1; j <= sel.GetNP(); j++)
		  {
			ppoints.Set(sel.PNum(j));
		  }
		  }
	  }

	  for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	  {
		  Element el = mesh[ei];
		  if (el.GetIndex() == k)
		  {
		  int todel = 0;
		  for (j = 0; j < el.GetNP(); j++)
		  {
			if (ppoints.Test(el[j]))
			{
			  todel = 1;
			}
		  }

		  if (el.GetNP() != 4)
		  {
			todel = 0;
		  }

		  if (todel != 0)
		  {
			  mesh[ei].Delete();
			  // ei--;
		  }
		  }
	  }
	  }

	  mesh.Compress();
	  PrintMessage(3, "Elements after Remove: ", mesh.GetNE());
	}

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//extern const char *ngscript[];
	public static void ExportArray<T, int BASE = 0, TIND = int>(py.module m)
	{
//C++ TO C# CONVERTER TODO TASK: There is no C# equivalent to the C++ 'typeid' operator:
	  string name = "Array_" + typeid(T).name();
  py.class_<Array<T,BASE,TIND>>(m, name).def("__len__", (Array<T,BASE,TIND> self) =>
  {
	  return self.Size();
  }).def("__getitem__", netgen.GlobalMembers.FunctionPointer((Array<T,BASE,TIND> self, TIND i) =>
  {
							 if (i < BASE || i >= BASE + self.Size())
							 {
							   throw py.index_error();
							 }
							 return self[i];
						   }), py.return_value_policy.reference).def("__iter__", (Array<T,BASE,TIND> self) =>
						   {
	return py.make_iterator(self.begin(),self.end());
	  }, py.keep_alive < 0,1>());
	}

	public static void TranslateException(NgException ex)
	{
	  string err = "Netgen exception: " + ex.What();
	  PyErr_SetString(PyExc_RuntimeError, err);
	}

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//static Transformation<3> global_trafo(Vec<3> (0,0,0));

	public static DLL_HEADER void ExportNetgenMeshing(py.module m)
	{
	  py.register_exception<NgException>(m, "NgException");
	  m.attr("_netgen_executable_started") = py.cast(netgen.netgen_executable_started);
	  string script;
//C++ TO C# CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged:
	  string * hcp = ngscript;
	  while (*hcp)
	  {
		  script += *hcp++;
	  }

	  m.attr("_ngscript") = py.cast(script);

  m.def("_GetStatus", () =>
  {
		  MyStr s = new MyStr();
		  double percent;
		  GetStatus(s.functorMethod, percent);
		  return py.make_tuple(s.c_str(), percent);
  });
  m.def("_PushStatus", (string s) =>
  {
	  PushStatus(new MyStr(s));
  });
  m.def("_SetThreadPercentage", (double percent) =>
  {
	  SetThreadPercent(percent);
  });

	//C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:
	  py.class_<NgMPI_Comm> (m, "MPI_Comm").def_property_readonly("rank", NgMPI_Comm.Rank).def_property_readonly("size", NgMPI_Comm.Size).def("Barrier", NgMPI_Comm.Barrier)

	#if PARALLEL
	.def("WTime", (NgMPI_Comm c) =>
	{
		return MPI_Wtime();
	} //C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:)#else.def("WTime", (NgMPI_Comm c) =>
	{
		return -1.0;
	} //C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:)#endif.def("Sum", (NgMPI_Comm c, double x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_SUM, c);
	}).def("Min", (NgMPI_Comm c, double x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_MIN, c);
	}).def("Max", (NgMPI_Comm c, double x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_MAX, c);
	}).def("Sum", (NgMPI_Comm c, int x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_SUM, c);
	}).def("Min", (NgMPI_Comm c, int x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_MIN, c);
	}).def("Max", (NgMPI_Comm c, int x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_MAX, c);
	}).def("Sum", (NgMPI_Comm c, uint x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_SUM, c);
	}).def("Min", (NgMPI_Comm c, uint x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_MIN, c);
	}).def("Max", (NgMPI_Comm c, uint x) =>
	{
		return MyMPI_AllReduceNG(x, MPI_MAX, c);
	}).def("SubComm", (NgMPI_Comm c, List<int> proc_list) =>
	{
		Array<int> procs = new Array<int>(proc_list.Count);
		for (int i = 0; i < procs.Size(); i++)
		{
		  procs[i] = proc_list[i];
		}
		if (!procs.Contains(c.Rank()))
		{
		  throw new Exception("rank " + ToString(c.Rank()) + " not in subcomm");
		}
	MPI_Comm subcomm = MyMPI_SubCommunicator(c, procs);
	return new NgMPI_Comm(subcomm, true);
	  }, py.arg("procs"));
	  ;




  py.class_<NGDummyArgument>(m, "NGDummyArgument").def("__bool__", (NGDummyArgument self) =>
  {
	  return false;
  });

  py.class_<Point < 2>> (m, "Point2d").def(py.init<double,double>()).def("__str__", ngcore.GlobalMembers.ToString<Point < 2>>).def(py.self - py.self).def(py.self + Vec < 2>()).def(py.self - Vec < 2>()).def("__getitem__", (Point < 2> self.functorMethod, int index) =>
  {
	  return self.functorMethod[index];
  });

  py.class_<Point < 3>> (m, "Point3d").def(py.init<double,double,double>()).def("__str__", ngcore.GlobalMembers.ToString<Point < 3>>).def(py.self - py.self).def(py.self + Vec < 3>()).def(py.self - Vec < 3>()).def("__getitem__", (Point < 2> self.functorMethod, int index) =>
  {
	  return self.functorMethod[index];
  });

  m.def("Pnt", netgen.GlobalMembers.FunctionPointer((double x, double y, double z) =>
  {
	  return global_trafo(Point < 3>(x,y,z));
  }));
  m.def("Pnt", netgen.GlobalMembers.FunctionPointer((double x, double y) =>
  {
	  return Point < 2>(x,y);
  }));

	  /*
	    // duplicated functions ????
	  m.def ("Pnt", FunctionPointer
	           ([](double x, double y, double z) { return Point<3>(x,y,z); }));
	  m.def ("Pnt", FunctionPointer
	           ([](double x, double y) { return Point<2>(x,y); }));
	  */

  py.class_<Vec < 2>> (m, "Vec2d").def(py.init<double,double>()).def("__str__", ngcore.GlobalMembers.ToString<Vec < 3>>).def(py.self + py.self).def(py.self - py.self).def(-py.self).def(double() * py.self).def("Norm", Vec < 2>.Length).def("__getitem__", (Vec < 2> vec.functorMethod, int index) =>
  {
	  return vec.functorMethod[index];
  }).def("__len__", (ref Vec < 2>) =>
  {
	  return 2;
  });

  py.class_<Vec < 3>> (m, "Vec3d").def(py.init<double,double,double>()).def("__str__", ngcore.GlobalMembers.ToString<Vec < 3>>).def(py.self + py.self).def(py.self - py.self).def(-py.self).def(double() * py.self).def("Norm", Vec < 3>.Length).def("__getitem__", (Vec < 3> vec.functorMethod, int index) =>
  {
	  return vec.functorMethod[index];
  }).def("__len__", (ref Vec < 3>) =>
  {
	  return 3;
  });

  m.def("Vec", netgen.GlobalMembers.FunctionPointer((double x, double y, double z) =>
  {
	  return global_trafo(Vec < 3>(x,y,z));
  }));
  m.def("Vec", netgen.GlobalMembers.FunctionPointer((double x, double y) =>
  {
	  return Vec < 2>(x,y);
  }));

  py.class_<Transformation < 3>> (m, "Trafo").def(py.init<Vec < 3>>(), "a translation").def(py.init<Point < 3>,Vec < 3>,double>(), "a rotation given by point on axes, direction of axes, angle").def("__mul__", (Transformation < 3> a.functorMethod, Transformation < 3> b.functorMethod) =>
  {
			 Transformation < 3> res;
			 res.Combine(a.functorMethod,b.functorMethod);
			 return res;
  }).def("__call__", (Transformation < 3> trafo.functorMethod, Point < 3> p.functorMethod) =>
  {
			 return new netgen.Transformation(trafo.functorMethod(p.functorMethod));
		 });

  m.def("GetTransformation", () =>
  {
	  return global_trafo;
  });
  m.def("SetTransformation", (Transformation < 3> trafo.functorMethod) =>
  {
	  global_trafo = trafo.functorMethod;
  });
  m.def("SetTransformation", (int dir, double angle) =>
  {
		   if (dir > 0)
		   {
			 global_trafo.SetAxisRotation(dir, angle * DefineConstants.M_PI / 180);
		   }
		   else
		   {
			 global_trafo = Transformation < 3> (Vec < 3>(0,0,0));
		   }
  }, py.arg("dir") = (int)0, py.arg("angle") = (int)0);
  m.def("SetTransformation", (Point < 3> p0.functorMethod, Vec < 3> ex.functorMethod, Vec < 3> ey.functorMethod, Vec < 3> ez.functorMethod) =>
  {
			  Point < 3> pnts[4];
			  pnts[0] = p0.functorMethod;
			  pnts[1] = p0.functorMethod + ex.functorMethod;
			  pnts[2] = p0.functorMethod + ey.functorMethod;
			  pnts[3] = p0.functorMethod + ez.functorMethod;
			  global_trafo = Transformation < 3> (pnts);
  }, py.arg("p0"), py.arg("ex"), py.arg("ey"), py.arg("ez"));



  py.class_<PointIndex>(m, "PointId").def(py.init<int>()).def("__repr__", ngcore.GlobalMembers.ToString<PointIndex>).def("__str__", ngcore.GlobalMembers.ToString<PointIndex>).def_property_readonly("nr", PointIndex.operator int).def("__eq__", netgen.GlobalMembers.FunctionPointer((PointIndex self, PointIndex other) =>
  {
					  return (int)self == (int)other;
  })).def("__hash__", netgen.GlobalMembers.FunctionPointer((PointIndex self) =>
  {
					  return (int)self;
				  }));

  py.class_<ElementIndex>(m, "ElementId3D").def(py.init<int>()).def("__repr__", ngcore.GlobalMembers.ToString<ElementIndex>).def("__str__", ngcore.GlobalMembers.ToString<ElementIndex>).def_property_readonly("nr", ElementIndex.operator int).def("__eq__", netgen.GlobalMembers.FunctionPointer((ElementIndex self, ElementIndex other) =>
  {
					  return (int)self == (int)other;
  })).def("__hash__", netgen.GlobalMembers.FunctionPointer((ElementIndex self) =>
  {
					  return (int)self;
				  }));


  py.class_<SurfaceElementIndex>(m, "ElementId2D").def(py.init<int>()).def("__repr__", ngcore.GlobalMembers.ToString<SurfaceElementIndex>).def("__str__", ngcore.GlobalMembers.ToString<SurfaceElementIndex>).def_property_readonly("nr", SurfaceElementIndex.operator int).def("__eq__", netgen.GlobalMembers.FunctionPointer((SurfaceElementIndex self, SurfaceElementIndex other) =>
  {
					  return (int)self == (int)other;
  })).def("__hash__", netgen.GlobalMembers.FunctionPointer((SurfaceElementIndex self) =>
  {
					  return (int)self;
				  }));

  py.class_<SegmentIndex>(m, "ElementId1D").def(py.init<int>()).def("__repr__", ngcore.GlobalMembers.ToString<SegmentIndex>).def("__str__", ngcore.GlobalMembers.ToString<SegmentIndex>).def_property_readonly("nr", SegmentIndex.operator int).def("__eq__", netgen.GlobalMembers.FunctionPointer((SegmentIndex self, SegmentIndex other) =>
  {
					  return (int)self == (int)other;
  })).def("__hash__", netgen.GlobalMembers.FunctionPointer((SegmentIndex self) =>
  {
					  return (int)self;
				  }));



	  /*  
	  py::class_<Point<3>> ("Point")
	    .def(py::init<double,double,double>())
	    ;
	  */

  py.class_<MeshPoint >(m, "MeshPoint").def(py.init<Point < 3>>()).def("__str__", ngcore.GlobalMembers.ToString<MeshPoint>).def("__repr__", ngcore.GlobalMembers.ToString<MeshPoint>).def_property_readonly("p", netgen.GlobalMembers.FunctionPointer((MeshPoint self) =>
  {
										 py.list l = new py.list();
										 l.append(py.cast(self[0]));
										 l.append(py.cast(self[1]));
										 l.append(py.cast(self[2]));
										 return py.tuple(l);
  })).def("__getitem__", netgen.GlobalMembers.FunctionPointer((MeshPoint self, int index) =>
  {
	  if (index < 0 || index>2)
	  {
			  throw py.index_error();
	  }
	  return self[index];
	})).def("__setitem__", netgen.GlobalMembers.FunctionPointer((MeshPoint self, int index, double val) =>
	{
	  if (index < 0 || index>2)
	  {
			  throw py.index_error();
	  }
	  self(index) = val;
	}));

  py.class_<Element>(m, "Element3D").def(py.init((int index, List<PointIndex> vertices) =>
  {
					int np = vertices.Count;
					ELEMENT_TYPE et;
					switch (np)
					{
					  case 4:
						  et = TET;
						  break;
					  case 5:
						  et = PYRAMID;
						  break;
					  case 6:
						  et = PRISM;
						  break;
					  case 8:
						  et = HEX;
						  break;
					  case 10:
						  et = TET10;
						  break;
					  case 13:
						  et = PYRAMID13;
						  break;
					  case 15:
						  et = PRISM15;
						  break;
					  case 20:
						  et = HEX20;
						  break;
					  default:
						throw new Exception("no Element3D with " + ToString(np) + " points");
					}

					var newel = new Element(et);
					for (int i = 0; i < np; i++)
					{
					  (*newel)[i] = vertices[i];
					}
					newel.SetIndex(index);
					return new netgen.Element(newel);
  }), py.arg("index") = 1,py.arg("vertices"), "create volume element").def("__repr__", ngcore.GlobalMembers.ToString<Element>).def_property("index", Element.GetIndex, Element.SetIndex).def_property("curved", Element.IsCurved, Element.SetCurved).def_property_readonly("vertices", netgen.GlobalMembers.FunctionPointer((Element self) =>
  {
									 py.list li = new py.list();
									 for (int i = 0; i < self.GetNV(); i++)
									 {
									   li.append(py.cast(self[i]));
									 }
									 return new py.list(li);
								   })).def_property_readonly("points", netgen.GlobalMembers.FunctionPointer((Element self) =>
								   {
									 py.list li = new py.list();
									 for (int i = 0; i < self.GetNP(); i++)
									 {
									   li.append(py.cast(self[i]));
									 }
									 return new py.list(li);
								   }));

  py.class_<Element2d>(m, "Element2D").def(py.init((int index, py.list vertices) =>
  {
					 Element2d newel = null;
					 if (py.len(vertices) == 3)
					 {
						 newel = new Element2d(TRIG);
						 for (int i = 0; i < 3; i++)
						 {
						   newel[i] = new py.extract<PointIndex>(vertices[i])();
						 }
						 newel.SetIndex(index);
					 }
					 else if (py.len(vertices) == 4)
					 {
						 newel = new Element2d(QUAD);
						 for (int i = 0; i < 4; i++)
						 {
						   newel[i] = new py.extract<PointIndex>(vertices[i])();
						 }
						 newel.SetIndex(index);
					 }
					 else if (py.len(vertices) == 6)
					 {
						 newel = new Element2d(TRIG6);
						 for (int i = 0; i < 6; i++)
						 {
						   newel[i] = new py.extract<PointIndex>(vertices[i])();
						 }
						 newel.SetIndex(index);
					 }
					 else if (py.len(vertices) == 8)
					 {
						 newel = new Element2d(QUAD8);
						 for (int i = 0; i < 8; i++)
						 {
						   newel[i] = new py.extract<PointIndex>(vertices[i])();
						 }
						 newel.SetIndex(index);
					 }
					 else
					 {
					   throw NgException("Inconsistent number of vertices in Element2D");
					 }
					 return newel;
  }), py.arg("index") = 1,py.arg("vertices"), "create surface element").def_property("index", Element2d.GetIndex, Element2d.SetIndex).def_property("curved", Element2d.IsCurved, Element2d.SetCurved).def_property_readonly("vertices", netgen.GlobalMembers.FunctionPointer((Element2d self) =>
  {
									py.list li = new py.list();
									for (int i = 0; i < self.GetNV(); i++)
									{
									  li.append(py.cast(self[i]));
									}
									return new py.list(li);
								  })).def_property_readonly("points", netgen.GlobalMembers.FunctionPointer((Element2d self) =>
								  {
									 py.list li = new py.list();
									 for (int i = 0; i < self.GetNP(); i++)
									 {
									   li.append(py.cast(self[i]));
									 }
									 return new py.list(li);
								   }));

  py.class_<Segment>(m, "Element1D").def(py.init((py.list vertices, py.list surfaces, int index, int edgenr) =>
  {
					Segment newel = new Segment();
					for (int i = 0; i < 2; i++)
					{
					  newel[i] = new py.extract<PointIndex>(vertices[i])();
					}
					newel.si = index;
					newel.edgenr = edgenr;
					newel.epgeominfo[0].edgenr = edgenr;
					newel.epgeominfo[1].edgenr = edgenr;
					// needed for codim2 in 3d
					newel.edgenr = index;
					if (len(surfaces))
					{
						newel.surfnr1 = new py.extract<int>(surfaces[0])();
						newel.surfnr2 = new py.extract<int>(surfaces[1])();
					}
					return newel;
  }), py.arg("vertices"), py.arg("surfaces") = py.list(), py.arg("index") = 1, py.arg("edgenr") = 1, "create segment element").def("__repr__", ngcore.GlobalMembers.ToString<Segment>).def_property_readonly("vertices", netgen.GlobalMembers.FunctionPointer((Segment self) =>
  {
									 py.list li = new py.list();
									 for (int i = 0; i < 2; i++)
									 {
									   li.append(py.cast(self[i]));
									 }
									 return new py.list(li);
								   })).def_property_readonly("points", netgen.GlobalMembers.FunctionPointer((Segment self) =>
								   {
									 py.list li = new py.list();
									 for (int i = 0; i < self.GetNP(); i++)
									 {
									   li.append(py.cast(self[i]));
									 }
									 return new py.list(li);
								   })).def_property_readonly("surfaces", netgen.GlobalMembers.FunctionPointer((Segment self) =>
								   {
									 py.list li = new py.list();
									 li.append(py.cast(self.surfnr1));
									 li.append(py.cast(self.surfnr2));
									 return new py.list(li);
								   })).def_property_readonly("index", netgen.GlobalMembers.FunctionPointer((Segment self) =>
								   {
			return self.si;
		  })).def_property_readonly("edgenr", netgen.GlobalMembers.FunctionPointer((Segment self) =>
		  {
							   return self.edgenr;
							 }));


  py.class_<Element0d>(m, "Element0D").def(py.init((PointIndex vertex, int index) =>
  {
					Element0d instance = new Element0d();
					instance.pnum.CopyFrom(vertex);
					instance.index = index;
					return instance;
  }), py.arg("vertex"), py.arg("index") = 1, "create point element").def("__repr__", ngcore.GlobalMembers.ToString<Element0d>).def_property_readonly("vertices", netgen.GlobalMembers.FunctionPointer((Element0d self) =>
  {
									 py.list li = new py.list();
									 li.append(py.cast(self.pnum));
									 return new py.list(li);
								   }));





  py.class_<FaceDescriptor>(m, "FaceDescriptor").def(py.init<const FaceDescriptor&>()).def(py.init((int surfnr, int domin, int domout, int bc) =>
  {
					FaceDescriptor instance = new FaceDescriptor();
					instance.SetSurfNr(surfnr);
					instance.SetDomainIn(domin);
					instance.SetDomainOut(domout);
					instance.SetBCProperty(bc);
							 return instance;
  }), py.arg("surfnr") = 1, py.arg("domin") = 1, py.arg("domout") = py.int_(0), py.arg("bc") = py.int_(0), "create facedescriptor").def("__str__", ngcore.GlobalMembers.ToString<FaceDescriptor>).def("__repr__", ngcore.GlobalMembers.ToString<FaceDescriptor>).def_property("surfnr", FaceDescriptor.SurfNr, FaceDescriptor.SetSurfNr).def_property("domin", FaceDescriptor.DomainIn, FaceDescriptor.SetDomainIn).def_property("domout", FaceDescriptor.DomainOut, FaceDescriptor.SetDomainOut).def_property("bc", FaceDescriptor.BCProperty, FaceDescriptor.SetBCProperty).def_property("bcname", (FaceDescriptor self) =>
  {
					  return self.GetBCName();
				  }, (FaceDescriptor self, string name) =>
				  {
					  self.SetBCName(new string(name));
				  } // memleak).def("SetSurfaceColor", (FaceDescriptor self, py.list color) =>
				  {
			Vec3d c = new Vec3d();
			c.X() = new py.extract<double>(color[0])();
			c.Y() = new py.extract<double>(color[1])();
			c.Z() = new py.extract<double>(color[2])();
			self.SetSurfColour(new netgen.Vec3d(c));
		  });



	  ExportArray<Element,0,uint>(m);
	  ExportArray<Element2d,0,uint>(m);
	  ExportArray<Segment,0,uint>(m);
	  ExportArray<Element0d>(m);
	  ExportArray<MeshPoint,PointIndex.BASE,PointIndex>(m);
	  ExportArray<FaceDescriptor>(m);

	  py.implicitly_convertible< int, PointIndex>();

	  py.class_<NetgenGeometry, NetgenGeometry> (m, "NetgenGeometry", py.dynamic_attr());

  py.class_<Mesh,Mesh>(m, "Mesh").def(py.init( (int dim, NgMPI_Comm comm) =>
  {
					 var mesh = new Mesh();
			 mesh.SetCommunicator(comm);
					 mesh.SetDimension(dim);
					 SetGlobalMesh(mesh); // for visualization
					 mesh.SetGeometry(null);
					 return new netgen.Mesh(mesh);
  }), py.arg("dim") = 3, py.arg("comm") = NgMPI_Comm{}).def(ngcore.GlobalMembers.NGSPickle<Mesh>()).def_property_readonly("comm", (Mesh amesh) =>
  {
				   return new ngcore.NgMPI_Comm(amesh.GetCommunicator());
			   }, "MPI-communicator the Mesh lives in").def_property_readonly("_timestamp", Mesh.GetTimeStamp).def("Distribute", (Mesh self, NgMPI_Comm comm) =>
			   {
	self.SetCommunicator(comm);
	if (comm.Size() == 1)
	{
		return self;
	}
	// if(MyMPI_GetNTasks(comm)==2) throw NgException("Sorry, cannot handle communicators with NP=2!");
	// cout << " rank " << MyMPI_GetId(comm) << " of " << MyMPI_GetNTasks(comm) << " called Distribute " << endl;
	if (comm.Rank() == 0)
	{
		self.Distribute();
	}
	else
	{
		self.SendRecvMesh();
	}
	return self;
	  }, py.arg("comm")).def("Receive", (NgMPI_Comm comm) =>
	  {
		var mesh = new Mesh();
		mesh.SetCommunicator(comm);
		mesh.SendRecvMesh();
		return new netgen.Mesh(mesh);
	  }).def("Load", netgen.GlobalMembers.FunctionPointer((Mesh self, string filename) =>
	  {

		var comm = self.GetCommunicator();
		int id = comm.Rank();
		int ntasks = comm.Size();
		var mesh = self;

		{
		  ifstream infile = new ifstream(filename);
		  if (!infile.good())
		  {
		throw NgException("Error opening file " + filename);
		  }
		}

		if (filename.IndexOf(".vol") == -1)
		{
		if (ntasks > 1)
		{
		  throw NgException("Not sure what to do with this?? Does this work with MPI??");
		}
		mesh.SetCommunicator(new ngcore.NgMPI_Comm(comm));
		ReadFile(*mesh,filename);
		//mesh->SetGlobalH (mparam.maxh);
		//mesh->CalcLocalH();
		return;
		}

		istream infile;
		Array<char> buf = new Array<char>(); // for distributing geometry!
		int strs;

		if (id == 0)
		{

		  if (filename.Substring(filename.Length - 3, 3) == ".gz")
		  {
		infile = new igzstream(filename);
		  }
		  else
		  {
		infile = new ifstream(filename);
		  }
		  mesh.Load(infile);

		  // make string from rest of file (for geometry info!)
		  // (this might be empty, in which case we take the global ng_geometry)
		  stringstream geom_part = new stringstream();
		  geom_part << infile.rdbuf();
		  string geom_part_string = geom_part.str();
		  strs = geom_part_string.Length;
		  // buf = new char[strs];
		  buf.SetSize(strs);
		  memcpy(buf[0], geom_part_string, strs * sizeof(char));

		  infile = null;

		  if (ntasks > 1)
		  {

		  string weightsfilename = new string(new char[filename.Length]);
		  weightsfilename = filename;
		  weightsfilename = StringFunctions.ChangeCharacter(weightsfilename, weightsfilename.Length - 3, 'w');
		  weightsfilename = StringFunctions.ChangeCharacter(weightsfilename, weightsfilename.Length - 2, 'e');
		  weightsfilename = StringFunctions.ChangeCharacter(weightsfilename, weightsfilename.Length - 1, 'i');

		  ifstream weightsfile = new ifstream(weightsfilename);
		  weightsfilename = null;

		  if (!(weightsfile.good()))
		  {
			  // cout << "regular distribute" << endl;
			  mesh.Distribute();
		  }
		  else
		  {
			  string str = new string(new char[20]);
			  bool endfile = false;
			  int n;
			  int dummy;

			  Array<int> segment_weights = new Array<int>();
			  Array<int> surface_weights = new Array<int>();
			  Array<int> volume_weights = new Array<int>();

			  while (weightsfile.good() && !endfile)
			  {
			  weightsfile >> str;

			  if (string.Compare(str, "edgeweights") == 0)
			  {
				  weightsfile >> n;
				  segment_weights.SetSize(n);
				  for (int i = 0; i < n; i++)
				  {
				weightsfile >> dummy >> segment_weights[i];
				  }
			  }

			  if (string.Compare(str, "surfaceweights") == 0)
			  {
				  weightsfile >> n;
				  surface_weights.SetSize(n);
				  for (int i = 0; i < n; i++)
				  {
				weightsfile >> dummy >> surface_weights[i];
				  }
			  }

			  if (string.Compare(str, "volumeweights") == 0)
			  {
				  weightsfile >> n;
				  volume_weights.SetSize(n);
				  for (int i = 0; i < n; i++)
				  {
				weightsfile >> dummy >> volume_weights[i];
				  }
			  }

			  if (string.Compare(str, "endfile") == 0)
			  {
				endfile = true;
			  }
			  }

			  mesh.Distribute(volume_weights, surface_weights, segment_weights);
		  }
		  } // ntasks>1 end
		} // id==0 end
		else
		{
		  mesh.SendRecvMesh();
		}

		if (ntasks > 1)
		{
//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to references to variables:
//ORIGINAL LINE: auto & mesh = self;
//C++ TO C# CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in C#:
#if PARALLEL
		  /** Scatter the geometry-string (no dummy-implementation in mpi_interface) **/
		  int strs = buf.Size();
		  MyMPI_Bcast(strs, comm);
		  if (strs > 0)
		  {
		MyMPI_Bcast(buf, comm);
		  }
//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to references to variables:
//ORIGINAL LINE: auto & mesh = self;
//C++ TO C# CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in C#:
#endif
		}

		NetgenGeometry geo;
		if (buf.Size())
		{ // if we had geom-info in the file, take it
		  istringstream geom_infile = new istringstream((string)(buf.Size(), (string) buf[0]));
		  geo = geometryregister.LoadFromMeshFile(geom_infile);
		}
		if (geo != null)
		{
			mesh.SetGeometry(geo);
		}
		else if (ng_geometry != null)
		{
			mesh.SetGeometry(ng_geometry);
		}
	  }),py.call_guard<py.gil_scoped_release>()).def("Save", static_cast < void(Mesh.)(string & name)const>(Mesh.Save), py.call_guard<py.gil_scoped_release>()).def("Export", (Mesh self, string filename, string format) =>
	  {
			if (WriteUserFormat(format, self, filename))
			{
				string err = "nothing known about format" + format;
				Array<char> names = new Array<char>();
				Array<char> extensions = new Array<char>();
				RegisterUserFormats(names, extensions);
				err += "\navailable formats are:\n";
				foreach (var name in names)
				{
				  err += "'" + name + "'\n";
				}
				throw NgException(err);
			}
		  }, py.arg("filename"), py.arg("format"),py.call_guard<py.gil_scoped_release>()).def_property("dim", Mesh.GetDimension, Mesh.SetDimension).def("Elements3D", (Array<Element,0,uint>&(Mesh.)()) (Mesh.VolumeElements), py.return_value_policy.reference).def("Elements2D", (Array<Element2d,0,uint>&(Mesh.)()) (Mesh.SurfaceElements), py.return_value_policy.reference).def("Elements1D", (Array<Segment,0,uint>&(Mesh.)()) (Mesh.LineSegments), py.return_value_policy.reference).def("Elements0D", netgen.GlobalMembers.FunctionPointer((Mesh self) =>
		  {
										 return new Array<Element0d>(self.pointelements);
									   }), py.return_value_policy.reference).def("Points", static_cast < Mesh.T_POINTS & (Mesh.)()> (Mesh.Points), py.return_value_policy.reference).def("FaceDescriptor", static_cast < FaceDescriptor & (Mesh.)(int)> (Mesh.GetFaceDescriptor), py.return_value_policy.reference).def("GetNFaceDescriptors", Mesh.GetNFD).def("GetNCD2Names", Mesh.GetNCD2Names).def("__getitem__", netgen.GlobalMembers.FunctionPointer((Mesh self, PointIndex pi) =>
									   {
										   return new netgen.MeshPoint(self[pi]);
										 })).def("Add", netgen.GlobalMembers.FunctionPointer((Mesh self, MeshPoint p) =>
										 {
									return new netgen.PointIndex(self.AddPoint(new Point3d(p)));
								  })).def("Add", netgen.GlobalMembers.FunctionPointer((Mesh self, Element el) =>
								  {
									return new netgen.ElementIndex(self.AddVolumeElement(el));
								  })).def("Add", netgen.GlobalMembers.FunctionPointer((Mesh self, Element2d el) =>
								  {
									return new netgen.SurfaceElementIndex(self.AddSurfaceElement(el));
								  })).def("Add", netgen.GlobalMembers.FunctionPointer((Mesh self, Segment el) =>
								  {
									return new netgen.SegmentIndex(self.AddSegment(el));
								  })).def("Add", netgen.GlobalMembers.FunctionPointer((Mesh self, Element0d el) =>
								  {
									return self.pointelements.Append(el);
								  })).def("Add", netgen.GlobalMembers.FunctionPointer((Mesh self, FaceDescriptor fd) =>
								  {
									return self.AddFaceDescriptor(fd);
								  })).def("DeleteSurfaceElement", netgen.GlobalMembers.FunctionPointer((Mesh self, SurfaceElementIndex i) =>
								  {
							 return self.DeleteSurfaceElement(i);
						   })).def("Compress", netgen.GlobalMembers.FunctionPointer((Mesh self) =>
						   {
										 return self.Compress();
									   }),py.call_guard<py.gil_scoped_release>()).def("SetBCName", Mesh.SetBCName).def("GetBCName", netgen.GlobalMembers.FunctionPointer((Mesh self, int bc) =>
									   {
										   return self.GetBCName(bc);
									   })).def("SetMaterial", Mesh.SetMaterial).def("GetMaterial", netgen.GlobalMembers.FunctionPointer((Mesh self, int domnr) =>
									   {
											 return new string(self.GetMaterial(domnr));
										 })).def("GetCD2Name", Mesh.GetCD2Name).def("SetCD2Name", Mesh.SetCD2Name).def("GetCD3Name", Mesh.GetCD3Name).def("SetCD3Name", Mesh.SetCD3Name).def("AddPointIdentification", (Mesh self, py.@object pindex1, py.@object pindex2, int identnr, int type) =>
										 {
				 if (new py.extract<PointIndex>(pindex1).check() && new py.extract<PointIndex>(pindex2).check() != null)
				 {
				 self.GetIdentifications().Add(new py.extract<PointIndex>(pindex1)(), new py.extract<PointIndex>(pindex2)(), identnr);
				 self.GetIdentifications().SetType(identnr, Identifications.ID_TYPE(type)); // type = 2 ... periodic
				 }
						   }, py.arg("pid1"), py.arg("pid2"), py.arg("identnr"), py.arg("type")).def("CalcLocalH", Mesh.CalcLocalH).def("SetMaxHDomain", (Mesh self, py.list maxhlist) =>
						   {
			Array<double> maxh = new Array<double>();
			foreach (var el in maxhlist)
			{
			  maxh.Append(py.cast<double>(el));
			}
			self.SetMaxHDomain(maxh);
		  }).def("GenerateVolumeMesh", (Mesh self, py.@object pymp) =>
		  {
			 Console.Write("generate vol mesh");
			 Console.Write("\n");

			 MeshingParameters mp = new MeshingParameters();
			 {
			   py.gil_scoped_acquire acquire = new py.gil_scoped_acquire();
			 if (new py.extract<MeshingParameters>(pymp).check() != null)
			 {
			   mp = new py.extract<MeshingParameters>(pymp)();
			 }
			 else
			 {
				 mp.optsteps3d = 5;
			 }
			 }
			 MeshVolume(mp, self);
			 OptimizeVolume(mp, self);
		   }, py.arg("mp") = new NGDummyArgument(),py.call_guard<py.gil_scoped_release>()).def("OptimizeVolumeMesh", netgen.GlobalMembers.FunctionPointer((Mesh self) =>
		   {
			MeshingParameters mp = new MeshingParameters();
			mp.optsteps3d = 5;
			OptimizeVolume(mp, self);
		  }),py.call_guard<py.gil_scoped_release>()).def("Refine", netgen.GlobalMembers.FunctionPointer((Mesh self) =>
		  {
			 if (self.GetGeometry() != null)
			 {
			   self.GetGeometry().GetRefinement().Refine(self);
			 }
			 else
			 {
			   Refinement().Refine(self);
			 }
			 self.UpdateTopology();
		   }),py.call_guard<py.gil_scoped_release>()).def("SecondOrder", netgen.GlobalMembers.FunctionPointer((Mesh self) =>
		   {
			 if (self.GetGeometry() != null)
			 {
			   self.GetGeometry().GetRefinement().MakeSecondOrder(self);
			 }
			 else
			 {
			   Refinement().MakeSecondOrder(self);
			 }
		   })).def("GetGeometry", (Mesh self) =>
		   {
			   return new netgen.NetgenGeometry(self.GetGeometry());
		   }).def("SetGeometry", (Mesh self, NetgenGeometry geo) =>
		   {
			 self.SetGeometry(geo);
		   }).def("BuildSearchTree", Mesh.BuildElementSearchTree, py.call_guard<py.gil_scoped_release>()).def("BoundaryLayer", netgen.GlobalMembers.FunctionPointer((Mesh self, int bc, py.list thicknesses, int volnr, py.list materials) =>
		   {
			 int n = py.len(thicknesses);
			 BoundaryLayerParameters blp = new BoundaryLayerParameters();

			 for (int i = 1; i <= self.GetNFD(); i++)
			 {
			   if (self.GetFaceDescriptor(i).BCProperty() == bc)
			   {
				   blp.surfid.Append(i);
			   }
			 }

			 Console.Write("add layer at surfaces: ");
			 Console.Write(blp.surfid);
			 Console.Write("\n");

			 blp.prismlayers = n;
			 blp.growthfactor = 1.0;

			 // find max domain nr
			 int maxind = 0;
			 for (ElementIndex ei = 0; ei < self.GetNE(); ei++)
			 {
			   maxind = Math.Max(maxind, self[ei].GetIndex());
			 }
			 Console.Write("maxind = ");
			 Console.Write(maxind);
			 Console.Write("\n");
			 for (int i = 0; i < n; i++)
			 {
				 blp.heights.Append(new py.extract<double>(thicknesses[i])());
				 blp.new_matnrs.Append(maxind + 1 + i);
				 self.SetMaterial(maxind + 1 + i, new py.extract<string>(materials[i])().c_str());
			 }
			 blp.bulk_matnr = volnr;
			 GenerateBoundaryLayer(self, blp);
		   })).def("BoundaryLayer", netgen.GlobalMembers.FunctionPointer((Mesh self, int bc, double thickness, int volnr, string material) =>
		   {
			 BoundaryLayerParameters blp = new BoundaryLayerParameters();

			 for (int i = 1; i <= self.GetNFD(); i++)
			 {
			   if (self.GetFaceDescriptor(i).BCProperty() == bc)
			   {
				   blp.surfid.Append(i);
			   }
			 }

			 Console.Write("add layer at surfaces: ");
			 Console.Write(blp.surfid);
			 Console.Write("\n");

			 blp.prismlayers = 1;
			 blp.hfirst = thickness;
			 blp.growthfactor = 1.0;

			 // find max domain nr
			 int maxind = 0;
			 for (ElementIndex ei = 0; ei < self.GetNE(); ei++)
			 {
			   maxind = Math.Max(maxind, self[ei].GetIndex());
			 }
			 Console.Write("maxind = ");
			 Console.Write(maxind);
			 Console.Write("\n");
			 self.SetMaterial(maxind + 1, material);
			 blp.new_matnr = maxind + 1;
			 blp.bulk_matnr = volnr;
			 GenerateBoundaryLayer(self, blp);
		   })).def("EnableTable", (Mesh self, string name, bool set) =>
		   {
			if (name == "edges")
			{
			  const_cast<MeshTopology&>(self.GetTopology()).SetBuildEdges(set);
			}
			if (name == "faces")
			{
			  const_cast<MeshTopology&>(self.GetTopology()).SetBuildFaces(set);
			}
		  }, py.arg("name"), py.arg("set") = true).def("Scale", netgen.GlobalMembers.FunctionPointer((Mesh self, double factor) =>
		  {
					 for (var i = 0; i < self.GetNP(); i++)
					 {
					   self.Point(i).Scale(factor);
					 }
				   }));

  m.def("ImportMesh", (string filename) =>
  {
						var mesh = new Mesh();
						ReadFile(*mesh, filename);
						return new netgen.Mesh(mesh);
  }, py.arg("filename"), @"Import mesh from other file format, supported file formats are:
 Neutral format (*.mesh, *.emt)
 Surface file (*.surf)
 Universal format (*.unv)
 Olaf format (*.emt)
 Tet format (*.tet)
 Pro/ENGINEER format (*.fnf)
");
	  py.enum_<MESHING_STEP>(m,"MeshingStep").value("MESHEDGES",MESHCONST_MESHEDGES).value("MESHSURFACE",MESHCONST_OPTSURFACE).value("MESHVOLUME",MESHCONST_OPTVOLUME);

  py.class_<MeshingParameters> (m, "MeshingParameters").def(py.init<>()).def(py.init((double maxh, bool quad_dominated, int optsteps2d, int optsteps3d, MESHING_STEP perfstepsend, int only3D_domain, string meshsizefilename, double grading, double curvaturesafety, double segmentsperedge) =>
  {
					MeshingParameters instance = new MeshingParameters();
					instance.maxh = maxh;
					instance.quad = (int)quad_dominated;
					instance.optsteps2d = optsteps2d;
					instance.optsteps3d = optsteps3d;
					instance.only3D_domain_nr = only3D_domain;
					instance.perfstepsend = (int)perfstepsend;
					instance.meshsizefilename = meshsizefilename;

					instance.grading = grading;
					instance.curvaturesafety = curvaturesafety;
					instance.segmentsperedge = segmentsperedge;
					return instance;
  }), py.arg("maxh") = 1000, py.arg("quad_dominated") = false, py.arg("optsteps2d") = 3, py.arg("optsteps3d") = 3, py.arg("perfstepsend") = MESHCONST_OPTVOLUME, py.arg("only3D_domain") = 0, py.arg("meshsizefilename") = "", py.arg("grading") = 0.3, py.arg("curvaturesafety") = 2, py.arg("segmentsperedge") = 1, "create meshing parameters").def("__str__", ngcore.GlobalMembers.ToString<MeshingParameters>).def_property("maxh", netgen.GlobalMembers.FunctionPointer((MeshingParameters mp) =>
  {
					  return mp.maxh;
				  }), netgen.GlobalMembers.FunctionPointer((MeshingParameters mp, double maxh) =>
				  {
					  return mp.maxh = maxh;
				  })).def("RestrictH", netgen.GlobalMembers.FunctionPointer((MeshingParameters mp, double x, double y, double z, double h) =>
				  {
			mp.meshsize_points.Append(new MeshingParameters.MeshSizePoint(Point < 3> (x,y,z), h));
		  }), py.arg("x"), py.arg("y"), py.arg("z"), py.arg("h"));

  m.def("SetTestoutFile", netgen.GlobalMembers.FunctionPointer((string filename) =>
  {
											   testout = null;
											   testout = new ofstream(filename);
  }));

  m.def("SetMessageImportance", netgen.GlobalMembers.FunctionPointer((int importance) =>
  {
													 int old = printmessage_importance;
													 printmessage_importance = importance;
													 return old;
  }));


	}

//C++ TO C# CONVERTER WARNING: The following constructor is declared outside of its associated class:
	public static PYBIND11_MODULE(libmesh UnnamedParameter, m UnnamedParameter2)
	{
	  ExportNetgenMeshing(m);
	}
	#endif








	/** Draws 2D rules.
	    Visual testing of 2D meshing rules */
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void DrawRules();

// A special function for Hermann Landes, Erlangen


	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define DLL_HEADER __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define __assume(cond)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_EXPORT __declspec(dllexport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API_IMPORT __declspec(dllimport)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_EXPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NGCORE_API NGCORE_API_IMPORT
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_INLINE inline
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max) { if ((value)<(min) || (value)>=(max)) throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max)); }
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_START(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_USER_END(n)
	//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
	//ORIGINAL LINE: #define VT_TRACER(n)

	/*
	
	  Very special implementations ..
	  
	 */


	///
	public static void CutOffAndCombine(Mesh mesh, Mesh othermesh)
	{
	  int i;
	  int j;
	  int nse = othermesh.GetNSE();
	  int onp = othermesh.GetNP();

	  int ne = mesh.GetNE();

	  PrintMessage(1, "other mesh has ", othermesh.GetNP(), " points, ", othermesh.GetNSE(), " surface elements.");

	  Array<Box3d> otherbounds = new Array<Box3d>(nse);
	  Box3d otherbox = new Box3d();

	  double maxh = 0;
	  for (i = 1; i <= nse; i++)
	  {
		  Element2d sel = othermesh.SurfaceElement(i);
		  sel.GetBox(othermesh.Points(), otherbounds.Elem(i));

		  double loch = othermesh.GetH(new othermesh.Point(sel.PNum(1)));
		  otherbounds.Elem(i).Increase(loch);
		  if (loch > maxh)
		  {
			  maxh = loch;
		  }
	  }

	  otherbox.SetPoint(new othermesh.Point(1));
	  for (i = 1; i <= othermesh.GetNP(); i++)
	  {
		otherbox.AddPoint(new othermesh.Point(i));
	  }
	  otherbox.Increase(maxh);

	  for (i = 1; i <= ne; i++)
	  {
		  Box3d box = new Box3d();
		  int remove = 0;

		  Element el = mesh.VolumeElement(i);
		  el.GetBox(mesh.Points(), box);

		  if (i % 10000 == 0)
		  {
//C++ TO C# CONVERTER TODO TASK: The cout 'flush' manipulator is not converted by C++ to C# Converter:
//ORIGINAL LINE: cout << "+" << flush;
		Console.Write("+");
		  }

		  if (box.Intersect(otherbox))
		  {
		  for (j = 1; j <= nse && !remove; j++)
		  {
			if (box.Intersect(otherbounds.Get(j)))
			{
			  remove = 1;
			}
		  }
		  }

		  if (remove != 0)
		  {
		mesh.VolumeElement(i).Delete();
		  }
	  }
	  Console.Write("\n");

	  BitArray connected = new BitArray(mesh.GetNP());
	  connected.Clear();
	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
		  Element2d el = mesh.SurfaceElement(i);
		  for (j = 1; j <= 3; j++)
		  {
		connected.Set(el.PNum(j));
		  }
	  }

	  bool changed;
	  do
	  {
		  changed = false;
		  for (i = 1; i <= mesh.GetNE(); i++)
		  {
		  Element el = mesh.VolumeElement(i);
		  int has = 0;
		  int hasnot = 0;
		  if (el[0] != null)
		  {
			  for (j = 0; j < 4; j++)
			  {
			  if (connected.Test(el[j]))
			  {
				has = 1;
			  }
			  else
			  {
				hasnot = 1;
			  }
			  }
			  if (has != 0 && hasnot != 0)
			  {
			  changed = true;
			  for (j = 0; j < 4; j++)
			  {
				connected.Set(el[j]);
			  }
			  }
		  }
		  }
//C++ TO C# CONVERTER TODO TASK: The cout 'flush' manipulator is not converted by C++ to C# Converter:
//ORIGINAL LINE: cout << "." << flush;
		  Console.Write(".");
	  } while (changed);
	  Console.Write("\n");

	  for (i = 1; i <= mesh.GetNE(); i++)
	  {
		  Element el = mesh.VolumeElement(i);
		  int hasnot = 0;
		  if (el[0] != null)
		  {
		  for (j = 0; j < 4; j++)
		  {
			  if (!connected.Test(el[j]))
			  {
			hasnot = 1;
			  }
		  }
		  if (hasnot != 0)
		  {
			mesh.VolumeElement(i).Delete();
		  }
		  }
	  }

	  mesh.Compress();

	  mesh.FindOpenElements();
	  BitArray locked = new BitArray(mesh.GetNP());
	  locked.Set();
	  for (i = 1; i <= mesh.GetNOpenElements(); i++)
	  {
		for (j = 1; j <= 3; j++)
		{
		  locked.Clear(mesh.OpenElement(i).PNum(j));
		}
	  }

	  for (PointIndex i = 1; i <= locked.Size(); i++)
	  {
		if (locked.Test(i))
		{
		mesh.AddLockedPoint(i);
		}
	  }




	  Array<PointIndex> pmat = new Array<PointIndex>(onp);

	  for (i = 1; i <= onp; i++)
	  {
		pmat.Elem(i) = mesh.AddPoint(new othermesh.Point(i));
	  }

	  int fnum = mesh.AddFaceDescriptor(new FaceDescriptor(0, 0, 1, 0));

	  for (i = 1; i <= othermesh.GetNSE(); i++)
	  {
		  Element2d tri = othermesh.SurfaceElement(i);
		  for (j = 1; j <= 3; j++)
		  {
		tri.PNum(j) = pmat.Get(tri.PNum(j));
		  }
		  tri.SetIndex(fnum);
		  mesh.AddSurfaceElement(tri);
	  }

	  for (i = 1; i <= onp; i++)
	  {
		mesh.AddLockedPoint(pmat.Elem(i));
	  }

	  mesh.CalcSurfacesOfNode();
	  mesh.CalcLocalH(0.3);
	}

	public static void HelmholtzMesh(Mesh mesh)
	{
	  int i;
	  double ri;
	  double ra;
	  double rinf;

	  Console.Write("ri = ");
	  ri = double.Parse(ConsoleInput.ReadToWhiteSpace(true));
	  Console.Write("ra = ");
	  ra = double.Parse(ConsoleInput.ReadToWhiteSpace(true));
	  Console.Write("rinf = ");
	  rinf = double.Parse(ConsoleInput.ReadToWhiteSpace(true));

	  double det = ri * ra * rinf - ri * ri * rinf;
	  double a = (ri - rinf) / det;
	  double b = (ri * ri - ra * rinf) / det;
	  for (i = 1; i <= mesh.GetNP(); i++)
	  {
		  Point < 3> & p = new mesh.Point(i);
		  double rold = ngsimd.GlobalMembers.sqrt(netgen.GlobalMembers.sqr(p(0)) + netgen.GlobalMembers.sqr(p(1)) + netgen.GlobalMembers.sqr(p(2)));
		  if (rold < ri)
		  {
			  continue;
		  }

		  double rnew = 1 / (a * rold - b);
		  double fac = rnew / rold;
		  p(0) *= fac;
		  p(1) *= fac;
		  p(2) *= fac;
	  }
	}
	public static HPREF_ELEMENT_TYPE ClassifyTet(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint)
	{
	  int ep1 = 0;
	  int ep2 = 0;
	  int ep3 = 0;
	  int ep4 = 0;
	  int cp1 = 0;
	  int cp2 = 0;
	  int cp3 = 0;
	  int cp4 = 0;
	  int fp1;
	  int fp2;
	  int fp3;
	  int fp4;
	  int isedge1 = 0;
	  int isedge2 = 0;
	  int isedge3 = 0;
	  int isedge4 = 0;
	  int isedge5 = 0;
	  int isedge6 = 0;
	  int isfedge1;
	  int isfedge2;
	  int isfedge3;
	  int isfedge4;
	  int isfedge5;
	  int isfedge6;
	  int isface1 = 0;
	  int isface2 = 0;
	  int isface3 = 0;
	  int isface4 = 0;

	  HPREF_ELEMENT_TYPE type = HP_NONE;


	  int debug = 0;
	  for (int j = 0;j < 4; j++)
	  {
		  if (el.pnums[j] == 444)
		  {
			  debug++;
		  }
		  if (el.pnums[j] == 115)
		  {
			  debug++;
		  }
		  if (el.pnums[j] == 382)
		  {
			  debug++;
		  }
		  if (el.pnums[j] == 281)
		  {
			  debug++;
		  }
	  }
	  if (debug < 4)
	  {
		  debug = 0;
	  }



	  for (int j = 0; j < 4; j++)
	  {
		for (int k = 0; k < 4; k++)
		{
		if (j == k)
		{
			continue;
		}
		if (type != null)
		{
			break;
		}

		int pi3 = 0;
		while (pi3 == j || pi3 == k)
		{
			pi3++;
		}
		int pi4 = 6 - j - k - pi3;

		// preserve orientation
		int[] sort = new int[4];
		sort[0] = j;
		sort[1] = k;
		sort[2] = pi3;
		sort[3] = pi4;
		int cnt = 0;
		for (int jj = 0; jj < 4; jj++)
		{
		  for (int kk = 0; kk < 3; kk++)
		  {
			if (sort[kk] > sort[kk + 1])
			{
			cnt++;
			Swap(sort[kk], sort[kk + 1]);
			}
		  }
		}
		if (cnt % 2 == 1)
		{
			Swap(pi3, pi4);
		}

		ep1 = edgepoint.Test(el.pnums[j]);
		ep2 = edgepoint.Test(el.pnums[k]);
		ep3 = edgepoint.Test(el.pnums[pi3]);
		ep4 = edgepoint.Test(el.pnums[pi4]);

		cp1 = cornerpoint.Test(el.pnums[j]);
		cp2 = cornerpoint.Test(el.pnums[k]);
		cp3 = cornerpoint.Test(el.pnums[pi3]);
		cp4 = cornerpoint.Test(el.pnums[pi4]);

		isedge1 = edges.Used(INDEX_2.Sort(el.pnums[j], el.pnums[k]));
		isedge2 = edges.Used(INDEX_2.Sort(el.pnums[j], el.pnums[pi3]));
		isedge3 = edges.Used(INDEX_2.Sort(el.pnums[j], el.pnums[pi4]));
		isedge4 = edges.Used(INDEX_2.Sort(el.pnums[k], el.pnums[pi3]));
		isedge5 = edges.Used(INDEX_2.Sort(el.pnums[k], el.pnums[pi4]));
		isedge6 = edges.Used(INDEX_2.Sort(el.pnums[pi3], el.pnums[pi4]));

		if (debug != 0)
		{
			Console.Write("debug");
			Console.Write("\n");
			*testout << "debug" << "\n";
			*testout << "ep = " << ep1 << ep2 << ep3 << ep4 << "\n";
			*testout << "cp = " << cp1 << cp2 << cp3 << cp4 << "\n";
			*testout << "edge = " << isedge1 << isedge2 << isedge3 << isedge4 << isedge5 << isedge6 << "\n";
		}


		isface1 = isface2 = isface3 = isface4 = 0;
		for (int l = 0; l < 4; l++)
		{
			INDEX_3 i3 = new INDEX_3(0,0,0);
			switch (l)
			{
				  case 0:
					  i3.I1() = el.pnums[k];
					  i3.I1() = el.pnums[pi3];
					  i3.I1() = el.pnums[pi4];
					  break;
				  case 1:
					  i3.I1() = el.pnums[j];
					  i3.I1() = el.pnums[pi3];
					  i3.I1() = el.pnums[pi4];
					  break;
				  case 2:
					  i3.I1() = el.pnums[j];
					  i3.I1() = el.pnums[k];
					  i3.I1() = el.pnums[pi4];
					  break;
				  case 3:
					  i3.I1() = el.pnums[j];
					  i3.I1() = el.pnums[k];
					  i3.I1() = el.pnums[pi3];
					  break;
			}
			i3.Sort();
			if (faces.Used(i3))
			{
			int domnr = faces.Get(i3);
			if (domnr == -1 || domnr == el.GetIndex())
			{
				switch (l)
				{
				  case 0:
					  isface1 = 1;
					  break;
				  case 1:
					  isface2 = 1;
					  break;
				  case 2:
					  isface3 = 1;
					  break;
				  case 3:
					  isface4 = 1;
					  break;
				}
			}
			}
		}
		/*
		  isface1 = faces.Used (INDEX_3::Sort (el.pnums[k], el.pnums[pi3], el.pnums[pi4]));
		  isface2 = faces.Used (INDEX_3::Sort (el.pnums[j], el.pnums[pi3], el.pnums[pi4]));
		  isface3 = faces.Used (INDEX_3::Sort (el.pnums[j], el.pnums[k], el.pnums[pi4]));
		  isface4 = faces.Used (INDEX_3::Sort (el.pnums[j], el.pnums[k], el.pnums[pi3]));
		*/

		isfedge1 = isfedge2 = isfedge3 = isfedge4 = isfedge5 = isfedge6 = 0;
		for (int l = 0; l < 6; l++)
		{
			INDEX_2 i2 = new INDEX_2(0,0);
			switch (l)
			{
				  case 0:
					  i2.I1() = el.pnums[j];
					  i2.I2() = el[k];
					  break;
				  case 1:
					  i2.I1() = el.pnums[j];
					  i2.I2() = el.pnums[pi3];
					  break;
				  case 2:
					  i2.I1() = el.pnums[j];
					  i2.I2() = el.pnums[pi4];
					  break;
				  case 3:
					  i2.I1() = el.pnums[k];
					  i2.I2() = el.pnums[pi3];
					  break;
				  case 4:
					  i2.I1() = el.pnums[k];
					  i2.I2() = el.pnums[pi4];
					  break;
				  case 5:
					  i2.I1() = el.pnums[pi3];
					  i2.I2() = el.pnums[pi4];
					  break;
			}
			i2.Sort();
			if (face_edges.Used(i2))
			{
			int domnr = face_edges.Get(i2);
			if (domnr == -1 || domnr == el.GetIndex())
			{
				switch (l)
				{
				  case 0:
					  isfedge1 = 1;
					  break;
				  case 1:
					  isfedge2 = 1;
					  break;
				  case 2:
					  isfedge3 = 1;
					  break;
				  case 3:
					  isfedge4 = 1;
					  break;
				  case 4:
					  isfedge5 = 1;
					  break;
				  case 5:
					  isfedge6 = 1;
					  break;
				}
			}
			}
		}
		/*
		  isfedge1 = face_edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[k]));
		  isfedge2 = face_edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[pi3]));
		  isfedge3 = face_edges.Used (INDEX_2::Sort (el.pnums[j], el.pnums[pi4]));
		  isfedge4 = face_edges.Used (INDEX_2::Sort (el.pnums[k], el.pnums[pi3]));
		  isfedge5 = face_edges.Used (INDEX_2::Sort (el.pnums[k], el.pnums[pi4]));
		  isfedge6 = face_edges.Used (INDEX_2::Sort (el.pnums[pi3], el.pnums[pi4]));
		*/

		fp1 = fp2 = fp3 = fp4 = 0;
		for (int l = 0; l < 4; l++)
		{
			int pti = 0;
			switch (l)
			{
			  case 0:
				  pti = el.pnums[j];
				  break;
			  case 1:
				  pti = el.pnums[k];
				  break;
			  case 2:
				  pti = el.pnums[pi3];
				  break;
			  case 3:
				  pti = el.pnums[pi4];
				  break;
			}
			int domnr = facepoint[pti];
			if (domnr == -1 || domnr == el.GetIndex())
			{
			switch (l)
			{
			  case 0:
				  fp1 = 1;
				  break;
			  case 1:
				  fp2 = 1;
				  break;
			  case 2:
				  fp3 = 1;
				  break;
			  case 3:
				  fp4 = 1;
				  break;
			}
			}
		}

		/*
		  fp1 = facepoint[el.pnums[j]] != 0;
		  fp2 = facepoint[el.pnums[k]] != 0;
		  fp3 = facepoint[el.pnums[pi3]] != 0;
		  fp4 = facepoint[el.pnums[pi4]] != 0;
		*/


		switch (isface1 + isface2 + isface3 + isface4)
		{
		  case 0:
		  {
			  isedge1 |= isfedge1;
			  isedge2 |= isfedge2;
			  isedge3 |= isfedge3;
			  isedge4 |= isfedge4;
			  isedge5 |= isfedge5;
			  isedge6 |= isfedge6;

			  ep1 |= fp1;
			  ep2 |= fp2;
			  ep3 |= fp3;
			  ep4 |= fp4;

			  switch (isedge1 + isedge2 + isedge3 + isedge4 + isedge5 + isedge6)
			  {
			case 0:
			{
				if (ep1 == 0 && ep2 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET;
				}

				if (ep1 != 0 && ep2 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_0E_1V;
				}

				if (ep1 != 0 && ep2 != 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_0E_2V;
				}

				if (ep1 != 0 && ep2 != 0 && ep3 != 0 && ep4 == 0)
				{
				  type = HP_TET_0E_3V;
				}

				if (ep1 != 0 && ep2 != 0 && ep3 != 0 && ep4 != 0)
				{
				  type = HP_TET_0E_4V;
				}

				break;
			}

			case 1:
			{
				if (isedge1 == 0)
				{
					break;
				}

				if (cp1 == 0 && cp2 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1E_0V;
				}

				if (cp1 != 0 && cp2 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1E_1VA;
				}

				if (cp1 == 0 && cp2 == 0 && ep3 == 0 && ep4 != 0)
				{
				  type = HP_TET_1E_1VB;
				}

				if (cp1 != 0 && cp2 != 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1E_2VA;
				}

				if (cp1 != 0 && cp2 == 0 && ep3 != 0 && ep4 == 0)
				{
				  type = HP_TET_1E_2VB;
				}

				if (cp1 != 0 && cp2 == 0 && ep3 == 0 && ep4 != 0)
				{
				  type = HP_TET_1E_2VC;
				}

				if (cp1 == 0 && cp2 == 0 && ep3 != 0 && ep4 != 0)
				{
				  type = HP_TET_1E_2VD;
				}

				if (cp1 != 0 && cp2 != 0 && ep3 != 0 && ep4 == 0)
				{
				  type = HP_TET_1E_3VA;
				}

				if (cp1 != 0 && cp2 == 0 && ep3 != 0 && ep4 != 0)
				{
				  type = HP_TET_1E_3VB;
				}

				if (cp1 != 0 && cp2 != 0 && ep3 != 0 && ep4 != 0)
				{
				  type = HP_TET_1E_4V;
				}

				break;
			}
			case 2:
			{
				if (isedge1 != 0 && isedge2 != 0)
				{
				if (cp2 == 0 && cp3 == 0 && ep4 == 0)
				{
				  type = HP_TET_2EA_0V;
				}

				if (cp2 != 0 && cp3 == 0 && ep4 == 0)
				{
				  type = HP_TET_2EA_1VA;
				}
				if (cp2 == 0 && cp3 != 0 && ep4 == 0)
				{
				  type = HP_TET_2EA_1VB;
				}

				if (cp2 == 0 && cp3 == 0 && ep4 != 0)
				{
				  type = HP_TET_2EA_1VC;
				}

				if (cp2 != 0 && cp3 != 0 && ep4 == 0)
				{
				  type = HP_TET_2EA_2VA;
				}
				if (cp2 != 0 && cp3 == 0 && ep4 != 0)
				{
				  type = HP_TET_2EA_2VB;
				}
				if (cp2 == 0 && cp3 != 0 && ep4 != 0)
				{
				  type = HP_TET_2EA_2VC;
				}

				if (cp2 != 0 && cp3 != 0 && ep4 != 0)
				{
				  type = HP_TET_2EA_3V;
				}
				}
				if (isedge1 != 0 && isedge6 != 0)
				{
				if (cp1 == 0 && cp2 == 0 && cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_2EB_0V;
				}
				if (cp1 != 0 && cp2 == 0 && cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_2EB_1V;
				}
				if (cp1 != 0 && cp2 != 0 && cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_2EB_2VA;
				}
				if (cp1 != 0 && cp2 == 0 && cp3 != 0 && cp4 == 0)
				{
				  type = HP_TET_2EB_2VB;
				}
				if (cp1 != 0 && cp2 == 0 && cp3 == 0 && cp4 != 0)
				{
				  type = HP_TET_2EB_2VC;
				}
				if (cp1 != 0 && cp2 != 0 && cp3 != 0 && cp4 == 0)
				{
				  type = HP_TET_2EB_3V;
				}
				if (cp1 != 0 && cp2 != 0 && cp3 != 0 && cp4 != 0)
				{
				  type = HP_TET_2EB_4V;
				}
				}
				break;
			}
			case 3:
			{
				if (isedge1 != 0 && isedge2 != 0 && isedge3 != 0)
				{
				if (cp2 == 0 && cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_3EA_0V;
				}
				if (cp2 != 0 && cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_3EA_1V;
				}
				if (cp2 != 0 && cp3 != 0 && cp4 == 0)
				{
				  type = HP_TET_3EA_2V;
				}
				if (cp2 != 0 && cp3 != 0 && cp4 != 0)
				{
				  type = HP_TET_3EA_3V;
				}
				}
				if (isedge1 != 0 && isedge3 != 0 && isedge4 != 0)
				{
				if (cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_3EB_0V;
				}
				if (cp3 != 0 && cp4 == 0)
				{
							  type = HP_TET_3EB_1V;
				}
				if (cp3 != 0 && cp4 != 0)
				{
				  type = HP_TET_3EB_2V;
				}
				}
				if (isedge1 != 0 && isedge2 != 0 && isedge5 != 0)
				{
				if (cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_3EC_0V;
				}
				if (cp3 != 0 && cp4 == 0)
				{
				  type = HP_TET_3EC_1V;
				}
				if (cp3 != 0 && cp4 != 0)
				{
				  type = HP_TET_3EC_2V;
				}
				}
				break;
			}
			  }
			  break;
		  }



		  case 1: // one singular face
		  {
			  if (isface1 == 0)
			  {
				  break;
			  }

			  switch (isfedge1 + isfedge2 + isfedge3 + isedge4 + isedge5 + isedge6)
			  {
			case 0:
			{
				if (fp1 == 0 && ep2 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1F_0E_0V;
				}
				if (fp1 != 0 && ep2 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1F_0E_1VB;
				}
				if (fp1 == 0 && ep2 != 0 && !ep3 & !ep4)
				{
				  type = HP_TET_1F_0E_1VA;
				}
				break;
			}
			case 1:
			{
				if (isfedge1 != 0)
				{
				if (ep1 == 0 && ep3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1F_1EA_0V;
				}
				}
				if (isedge4 != 0) // V1-V3
				{
				if (ep1 == 0 && cp2 == 0 && cp3 == 0 && ep4 == 0)
				{
				  type = HP_TET_1F_1EB_0V;
				}
				}
				break;
			}
			  }
			  break;
		  }


		  case 2: // two singular faces
		  {
			  if (isface1 == 0 || isface2 == 0)
			  {
				  break;
			  }

			  switch (isfedge1 + isedge2 + isedge3 + isedge4 + isedge5)
			  {
			case 0:
			{
				if (ep1 == 0 && ep2 == 0 && cp3 == 0 && cp4 == 0)
				{
				  type = HP_TET_2F_0E_0V;
				}
				break;
			}
			  }
			  break;
		  }


		}

		if (type != HP_NONE)
		{
			int[] pnums = new int[4];
			pnums[0] = el.pnums[j];
			pnums[1] = el.pnums[k];
			pnums[2] = el.pnums[pi3];
			pnums[3] = el.pnums[pi4];
			for (k = 0;k < 4;k++)
			{
				el.pnums[k] = pnums[k];
			}
			break;
		}
		}
	  }


	  if (debug != 0)
	  {
		  Console.Write("type = ");
		  Console.Write(type);
		  Console.Write("\n");
	  }

	  if (type == HP_NONE)
	  {
		  //     cnt_undef++;
		  (*testout) << "undefined element" << "\n" << "cp = " << cp1 << cp2 << cp3 << cp4 << "\n" << "ep = " << ep1 << ep2 << ep3 << ep4 << "\n" << "isedge = " << isedge1 << isedge2 << isedge3 << isedge4 << isedge5 << isedge6 << "\n" << "isface = " << isface1 << isface2 << isface3 << isface4 << "\n";
		  Console.Write("undefined element !!! ");
		  Console.Write("\n");


	  }
	  return (type);
	}



	public static HPREF_ELEMENT_TYPE ClassifyPrism(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint)
	{

	  HPREF_ELEMENT_TYPE type = HP_NONE;

	  int[] p = new int[6];
	  for (int m = 1;m <= 6;m++)
	  {
		  int[] point_sing = {0, 0, 0, 0, 0, 0};
		  int[] face_sing = {0, 0, 0, 0, 0};
		  int[] edge_sing = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		  if (m < 4)
		  {
		  p[0] = m;
		  p[1] = m % 3 + 1;
		  p[2] = (m % 3 + 1) % 3 + 1;
		  for (int l = 3;l < 6;l++)
		  {
			  p[l] = p[l - 3] + 3;
		  }
		  }
		  else
		  {
		  p[0] = m;
		  p[1] = (m % 3 + 1) % 3 + 4;
		  p[2] = m % 3 + 4;
		  for (int l = 3;l < 6;l++)
		  {
			  p[l] = p[l - 3] - 3;
		  }
		  }

		  for (int j = 0;j < 6;j++)
		  {
		  if (cornerpoint.Test(el.PNum(p[j])))
		  {
			  point_sing[p[j] - 1] = 3;
		  }
		  else if (edgepoint.Test(el.PNum(p[j])))
		  {
			  point_sing[p[j] - 1] = 2;
		  }
		  else if (facepoint[el.PNum(p[j])] == -1 || facepoint[el.PNum(p[j])] == el.GetIndex())
		  {
			point_sing[p[j] - 1] = 1;
		  }
		  }

		  ELEMENT_EDGE[] eledges = MeshTopology.GetEdges1(PRISM);
		  for (int k = 0;k < 9;k++)
		  {
		  INDEX_2 i2 = INDEX_2.Sort(el.PNum(p[eledges[k][0] - 1]),el.PNum(p[eledges[k][1] - 1]));
		  if (edges.Used(i2))
		  {
			  edge_sing[k] = 2;
		  }
		  else
		  {
			  edge_sing[k] = face_edges.Used(i2);
		  }
		  }

		  ELEMENT_FACE[] elfaces = MeshTopology.GetFaces1(PRISM);
		  for (int k = 0;k < 5;k++)
		  {
		  INDEX_3 i3 = new INDEX_3();

		  if (k < 2)
		  {
			i3 = INDEX_3.Sort(el.pnums[p[elfaces[k][0] - 1] - 1], el.pnums[p[elfaces[k][1] - 1] - 1], el.pnums[p[elfaces[k][2] - 1] - 1]);
		  }
		  else
		  {
			  INDEX_4 i4 = new INDEX_4(el.pnums[p[elfaces[k][0] - 1] - 1], el.pnums[p[elfaces[k][1] - 1] - 1], el.pnums[p[elfaces[k][2] - 1] - 1],el.pnums[p[elfaces[k][3] - 1] - 1]);
			  i4.Sort();
			  i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3());
		  }

		  if (faces.Used(i3))
		  {
			  int domnr = faces.Get(i3);
			  if (domnr == -1 || domnr == el.GetIndex())
			  {
			face_sing[k] = 1;
			  }

		  }
		  }
		  if (face_sing[1] > face_sing[0])
		  {
			  m = m + 2;
			  continue;
		  }


		  //int cp = 0;  

		  int qfsing = face_sing[2] + face_sing[3] + face_sing[4];
		  int tfsing = face_sing[0] + face_sing[1];
		  int evsing = edge_sing[6] + edge_sing[7] + edge_sing[8];
		  int ehsing = edge_sing[0] + edge_sing[1] + edge_sing[2] + edge_sing[3] + edge_sing[4] + edge_sing[5];

		  if (qfsing + tfsing + evsing + ehsing == 0)
		  {
			type = HP_PRISM;
			break;
		  }

		  HPREF_ELEMENT_TYPE[] types = {HP_NONE, HP_NONE, HP_NONE};

		  int fb = (1 - face_sing[4]) * face_sing[3] * (face_sing[2] + face_sing[3]) + 3 * face_sing[4] * face_sing[3] * face_sing[2];
		  int[] sve = {edge_sing[7], edge_sing[8], edge_sing[6]};


		  if (fb != qfsing)
		  {
			  continue;
		  }


		  switch (fb)
		  {
		case 0:
		  if (evsing == 0 && ehsing == 3 * tfsing)
		  {
			  types[0] = HP_PRISM;
			  types[1] = HP_PRISM_1FA_0E_0V;
			  types[2] = HP_PRISM_2FA_0E_0V;
		  }
		  if (evsing > 0 && sve[0] == evsing) // 1 vertical edge 1-4
		  {
			  types[0] = HP_PRISM_SINGEDGE;
			  types[1] = HP_PRISM_1FA_1E_0V;
			  types[2] = HP_PRISM_2FA_1E_0V;
		  }

		  if (sve[0] > 0 && sve[1] > 0 && sve[2] == 0)
		  {
			  types[0] = HP_PRISM_SINGEDGE_V12;
			  types[1] = HP_PRISM_1FA_2E_0V;
			  types[2] = HP_PRISM_2FA_2E_0V;
		  }
		  if (sve[0] > 0 && sve[1] > 0 && sve[2] > 0)
		  {
			  types[0] = HP_PRISM_3E_0V;
			  types[1] = HP_PRISM_1FA_3E_0V;
			  types[2] = HP_PRISM_2FA_3E_0V;

			  if (edge_sing[0] > 1 && edge_sing[2] > 1 && edge_sing[4] > 1 && edge_sing[5] > 1 && tfsing == 0)
			  {
			types[0] = HP_PRISM_3E_4EH;
			  }
		  }

		  break;
		case 1:
		  if (sve[0] <= 1 && sve[1] <= 1)
		  {
				  if (sve[2] == 0)
				  {
					  types[0] = HP_PRISM_1FB_0E_0V;
					  types[1] = HP_PRISM_1FA_1FB_0E_0V;
					  types[2] = HP_PRISM_2FA_1FB_0E_0V;
				  }
				  else
				  {
					  types[0] = HP_PRISM_1FB_1EC_0V;
					  types[1] = HP_PRISM_1FA_1FB_1EC_0V;
					  types[2] = HP_PRISM_2FA_1FB_1EC_0V;
				  }
		  }

		  if (sve[0] > 1 && sve[2] >= 1 && sve[1] <= 1)
		  {
			  types[0] = HP_PRISM_1FB_2EB_0V;
			  types[1] = HP_PRISM_1FA_1FB_2EB_0V;
			  types[2] = HP_PRISM_2FA_1FB_2EB_0V;
		  }

		  if (sve[0] > 1 && sve[1] <= 1 && sve[2] == 0) // ea && !eb
		  {
			  types[0] = HP_PRISM_1FB_1EA_0V;
			  types[1] = HP_PRISM_1FA_1FB_1EA_0V;
			  types[2] = HP_PRISM_2FA_1FB_1EA_0V;
		  }

		  if (sve[0] <= 1 && sve[1] > 1 && sve[2] == 0)
		  {
			types[1] = HP_PRISM_1FA_1FB_1EB_0V;
		  }

		  if (sve[0] > 1 && sve[1] > 1)
		  {
			if (sve[2] == 0) // ea && eb
			{
			types[0] = HP_PRISM_1FB_2EA_0V;
			types[1] = HP_PRISM_1FA_1FB_2EA_0V;
			types[2] = HP_PRISM_2FA_1FB_2EA_0V;
			}
		  }
		  if (sve[0] <= 1 && sve[1] > 1 && sve[2] > 0)
		  {
			types[1] = HP_PRISM_1FA_1FB_2EC_0V;
		  }

		  if (sve[0] > 1 && sve[1] > 1 && sve[2] >= 1) //sve[2] can also be a face-edge
		  {
			  types[0] = HP_PRISM_1FB_3E_0V;
			  types[1] = HP_PRISM_1FA_1FB_3E_0V;
			  types[2] = HP_PRISM_2FA_1FB_3E_0V;
		  }

		  break;

		case 2:
		  if (sve[0] <= 1)
		  {
			Console.Write(" **** WARNING: Edge between to different singular faces should be marked singular ");
			Console.Write("\n");
		  }

		  if (sve[1] <= 1)
		  {
			if (sve[2] <= 1)
			{
			types[0] = HP_PRISM_2FB_0E_0V;
			types[1] = HP_PRISM_1FA_2FB_0E_0V;
			types[2] = HP_PRISM_2FA_2FB_0E_0V;
			}
			else
			{
			types[0] = HP_PRISM_2FB_1EC_0V;
			types[1] = HP_PRISM_1FA_2FB_1EC_0V;
			types[2] = HP_PRISM_2FA_2FB_1EC_0V;
			}
		  }
		  else
		  {
			if (sve[2] <= 1)
			{
			  types[1] = HP_PRISM_1FA_2FB_1EB_0V;
			}
			else
			{
			types[0] = HP_PRISM_2FB_3E_0V;
			types[1] = HP_PRISM_1FA_2FB_3E_0V;
			types[2] = HP_PRISM_2FA_2FB_3E_0V;
			}
		  }

		  break;

		case 3:
		  types[0] = HP_PRISM_3FB_0V;
		  types[1] = HP_PRISM_1FA_3FB_0V;
		  types[2] = HP_PRISM_2FA_3FB_0V;
		  break;
		  }
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: type = types[tfsing];
		  type.CopyFrom(types[tfsing]);


		  if (type != HP_NONE)
		  {
		break;
		  }
	  }

	  /*
	   *testout << " Prism with pnums " << endl; 
	   for(int j=0;j<6;j++) *testout << el.pnums[j] << "\t"; 
	   *testout << endl; 
	   */

	  if (type != HP_NONE)
	  {
		  int[] pnums = new int[6];
		  for (int j = 0;j < 6;j++)
		  {
			  pnums[j] = el.PNum(p[j]);
		  }
		  for (int k = 0;k < 6;k++)
		  {
			  el.pnums[k] = pnums[k];
		  }
	  }

	  /* *testout << " Classified Prism with pnums " << endl; 
	     for(int j=0;j<6;j++) *testout << el.pnums[j] << "\t"; 
	     *testout << endl; 
	     */ 
	  return (type);
	}


	// #ifdef SABINE 
	public static HPREF_ELEMENT_TYPE ClassifyTrig(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint, int dim, FaceDescriptor fd)

	{
	  HPREF_ELEMENT_TYPE type = HP_NONE;

	  int[] pnums = new int[3];
	  int[] p = new int[3];

	  INDEX_3 i3 = new INDEX_3(el.pnums[0], el.pnums[1], el.pnums[2]);
	  i3.Sort();
	  bool sing_face = faces.Used(i3);

	  // *testout << " facepoint " << facepoint << endl;  


	  // Try all rotations of the trig 
	  for (int j = 0;j < 3;j++)
	  {
		  int[] point_sing = {0, 0, 0};
		  int[] edge_sing = {0, 0, 0};
		  // *testout << " actual rotation of trig points " ;  
		  for (int m = 0;m < 3;m++)
		  {
		  p[m] = (j + m) % 3 + 1; // local vertex number
		  pnums[m] = el.PNum(p[m]); // global vertex number
		  // *testout << pnums[m] << " \t "; 
		  }
		  // *testout << endl ; 

		  if (dim == 3)
		  {
		  // face point 
		  for (int k = 0;k < 3;k++)
		  {
			if (!sing_face)
			{
			//	*testout << " fp [" << k << "] = " << facepoint[pnums[k]] << endl;   
			//	*testout << " fd.DomainIn()" <<  fd.DomainIn() << endl; 
			//	*testout  << " fd.DomainOut()" <<  fd.DomainOut() << endl; 
			if (facepoint[pnums[k]] && (facepoint[pnums[k]] == -1 || facepoint[pnums[k]] == fd.DomainIn() || facepoint[pnums[k]] == fd.DomainOut()))
			{
			  point_sing[p[k] - 1] = 1;
			}
			}
		  }
		  // if point is on face_edge in next step sing = 2 

		  /*	  *testout << " pointsing NACH FACEPOints ... FALLS EDGEPOINT UMSETZEN" ; 
				for (int k=0;k<3;k++) *testout << "\t" << point_sing[p[k]-1] ;
				*testout << endl; */
		  }

		  ELEMENT_EDGE[] eledges = MeshTopology.GetEdges1(TRIG);

		  if (dim == 3)
		  {
		  for (int k = 0;k < 3;k++)
		  {
			  int ep1 = p[eledges[k][0] - 1];
			  int ep2 = p[eledges[k][1] - 1];
			  INDEX_2 i2 = new INDEX_2(el.PNum(ep1),el.PNum(ep2));

			  if (edges.Used(i2))
			  {

			  edge_sing[k] = 2;
			  point_sing[ep1 - 1] = 2;
			  point_sing[ep2 - 1] = 2;
			  }
			  else // face_edge?
			  {
			  i2.Sort();
			  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1) // edge not face_edge acc. to surface in which trig lies
			  {
				  if (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut())
				  {
				  edge_sing[k] = 1;
				  }
				  else
				  {
				  point_sing[ep1 - 1] = 0; // set to edge_point
				  point_sing[ep2 - 1] = 0; // set to edge_point
				  }
			  }
			  }

			  /*  *testout << " pointsing NACH edges UND FACEEDGES UMSETZEN ... " ; 
					  for (int k=0;k<3;k++) *testout << "\t" << point_sing[p[k]-1] ; 
					  *testout << endl;          
					  */
		  }
		  }
		  /*
		   *testout << " dim " << dim << endl; 
		   *testout << " edgepoint_dom " << edgepoint_dom << endl; 
		   */
		  if (dim == 2)
		  {
		  for (int k = 0;k < 3;k++)
		  {
			  int ep1 = p[eledges[k][0] - 1];
			  int ep2 = p[eledges[k][1] - 1];

			  INDEX_2 i2 = new INDEX_2(el.PNum(ep1),el.PNum(ep2));

			  if (edges.Used(i2))
			  {

			  if (edgepoint_dom.Used(INDEX_2(fd.SurfNr(),pnums[ep1 - 1])) || edgepoint_dom.Used(INDEX_2(-1,pnums[ep1 - 1])) || edgepoint_dom.Used(INDEX_2(fd.SurfNr(),pnums[ep2 - 1])) || edgepoint_dom.Used(INDEX_2(-1,pnums[ep2 - 1])))
			  {
				  edge_sing[k] = 2;
				  point_sing[ep1 - 1] = 2;
				  point_sing[ep2 - 1] = 2;
			  }
			  }

		  }
		  }



		  for (int k = 0;k < 3;k++)
		  {
		if (edgepoint.Test(pnums[k])) //edgepoint, but not member of sing_edge on trig -> cp
		{
			INDEX_2 i2a = INDEX_2.Sort(el.PNum(p[k]), el.PNum(p[(k + 1) % 3]));
			INDEX_2 i2b = INDEX_2.Sort(el.PNum(p[k]), el.PNum(p[(k + 2) % 3]));

			if (!edges.Used(i2a) && !edges.Used(i2b))
			{
			  point_sing[p[k] - 1] = 3;
			}
		}
		  }

		  for (int k = 0;k < 3;k++)
		  {
		if (cornerpoint.Test(el.PNum(p[k])))
		{
		  point_sing[p[k] - 1] = 3;
		}
		  }

		  *testout << "point_sing = " << point_sing[0] << point_sing[1] << point_sing[2] << "\n";

		  if (edge_sing[0] + edge_sing[1] + edge_sing[2] == 0)
		  {
			  int ps = point_sing[0] + point_sing[1] + point_sing[2];

			  if (ps == 0)
			  {
				type = HP_TRIG;
			  }
			  else if (point_sing[p[0] - 1] && !point_sing[p[1] - 1] && !point_sing[p[2] - 1])
			  {
				type = HP_TRIG_SINGCORNER;
			  }
			  else if (point_sing[p[0] - 1] && point_sing[p[1] - 1] && !point_sing[p[2] - 1])
			  {
				type = HP_TRIG_SINGCORNER12;
			  }
			  else if (point_sing[p[0] - 1] && point_sing[p[1] - 1] && point_sing[p[2] - 1])
			  {
				  if (dim == 2)
				  {
					  type = HP_TRIG_SINGCORNER123_2D;
				  }
				  else
				  {
					  type = HP_TRIG_SINGCORNER123;
				  }
			  }
		  }
		  else
		  {
			if (edge_sing[2] != 0 && edge_sing[0] == 0 && edge_sing[1] == 0) //E[2]=(1,2)
			{
				int code = 0;
				if (point_sing[p[0] - 1] > edge_sing[2])
				{
					code += 1;
				}
				if (point_sing[p[1] - 1] > edge_sing[2])
				{
					code += 2;
				}
				if (point_sing[p[2] - 1] != 0)
				{
					code += 4;
				}

				HPREF_ELEMENT_TYPE[] types = {HP_TRIG_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER12, HP_TRIG_SINGEDGECORNER3, HP_TRIG_SINGEDGECORNER13, HP_TRIG_SINGEDGECORNER23, HP_TRIG_SINGEDGECORNER123};
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: type = types[code];
				type.CopyFrom(types[code]);

			} // E[0] = [0,2], E[1] =[1,2], E[2] = [0,1]
			else
			{
			  if (edge_sing[2] != 0 && edge_sing[1] == 0 && edge_sing[0] != 0)
			  {
				  if (point_sing[p[2] - 1] <= edge_sing[0])
				  {
					  if (point_sing[p[1] - 1] <= edge_sing[2])
					  {
						  type = HP_TRIG_SINGEDGES;
					  }
					  else
					  {
						  type = HP_TRIG_SINGEDGES2;
					  }
				  }
				  else
				  {
					  if (point_sing[p[1] - 1] <= edge_sing[2])
					  {
						type = HP_TRIG_SINGEDGES3;
					  }
					  else
					  {
						  type = HP_TRIG_SINGEDGES23;
					  }
				  }
			  }
			  else if (edge_sing[2] && edge_sing[1] && edge_sing[0])
			  {
				type = HP_TRIG_3SINGEDGES;
			  }
			}
		  }

		  //  cout << " run for " <<  j << " gives type " << type << endl; 
		  //*testout << " run for " <<  j << " gives type " << type << endl; 

		  if (type != HP_NONE)
		  {
			  break;
		  }
	  }

	  *testout << "type = " << type << "\n";

	  for (int k = 0;k < 3;k++)
	  {
		  el[k] = pnums[k];
	  }
	  /*if(type != HP_NONE) 
	    {
	     
	    cout << " TRIG with pnums " << pnums[0] << "\t"  << 
	    pnums[1] << "\t"  << pnums[2] << endl; 
	    cout << " type "  << type << endl; 
	    }
	  */
		  return (type);
	}
	#if HPREF_OLD
	public static HPREF_ELEMENT_TYPE ClassifyTrig(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint, int dim, FaceDescriptor fd)
	{
	  HPREF_ELEMENT_TYPE type = HP_NONE;

	  int[] pnums = new int[3];

	  INDEX_3 i3 = new INDEX_3(el.pnums[0], el.pnums[1], el.pnums[2]);
	  i3.Sort();
	  bool sing_face = faces.Used(i3);


	  for (int j = 1; j <= 3; j++)
	  {
		  int ep1 = edgepoint.Test(el.PNumMod(j));
		  int ep2 = edgepoint.Test(el.PNumMod(j + 1));
		  int ep3 = edgepoint.Test(el.PNumMod(j + 2));

		  if (dim == 2)
		  {
		  // JS, Dec 11
		  ep1 = edgepoint_dom.Used(INDEX_2(fd.SurfNr(), el.PNumMod(j))) || edgepoint_dom.Used(INDEX_2(-1, el.PNumMod(j)));
		  ep2 = edgepoint_dom.Used(INDEX_2(fd.SurfNr(), el.PNumMod(j + 1))) || edgepoint_dom.Used(INDEX_2(-1, el.PNumMod(j + 1)));
		  ep3 = edgepoint_dom.Used(INDEX_2(fd.SurfNr(), el.PNumMod(j + 2))) || edgepoint_dom.Used(INDEX_2(-1, el.PNumMod(j + 2)));
		  /*
				ep1 = edgepoint_dom.Used (INDEX_2 (el.index, el.PNumMod(j)));
				ep2 = edgepoint_dom.Used (INDEX_2 (el.index, el.PNumMod(j+1)));
				ep3 = edgepoint_dom.Used (INDEX_2 (el.index, el.PNumMod(j+2)));
		  */
		  // ep3 = edgepoint_dom.Used (INDEX_2 (mesh.SurfaceElement(i).GetIndex(), el.PNumMod(j+2)));
		  }



		  int cp1 = cornerpoint.Test(el.PNumMod(j));
		  int cp2 = cornerpoint.Test(el.PNumMod(j + 1));
		  int cp3 = cornerpoint.Test(el.PNumMod(j + 2));

		  ep1 |= cp1;
		  ep2 |= cp2;
		  ep3 |= cp3;


		  // (*testout) << "cp = " << cp1 << cp2 << cp3 << ", ep = " << ep1 << ep2 << ep3 << endl;

		  int[] p = {el.PNumMod(j), el.PNumMod(j + 1), el.PNumMod(j + 2)};
		  if (ep1 != 0)
		  {
		  INDEX_2 i2a = INDEX_2.Sort(p[0], p[1]);
		  INDEX_2 i2b = INDEX_2.Sort(p[0], p[2]);
		  if (!edges.Used(i2a) && !edges.Used(i2b))
		  {
			cp1 = 1;
		  }
		  }
		  if (ep2 != 0)
		  {
		  INDEX_2 i2a = INDEX_2.Sort(p[1], p[0]);
		  INDEX_2 i2b = INDEX_2.Sort(p[1], p[2]);
		  if (!edges.Used(i2a) && !edges.Used(i2b))
		  {
			cp2 = 1;
		  }
		  }
		  if (ep3 != 0)
		  {
		  INDEX_2 i2a = INDEX_2.Sort(p[2], p[0]);
		  INDEX_2 i2b = INDEX_2.Sort(p[2], p[1]);
		  if (!edges.Used(i2a) && !edges.Used(i2b))
		  {
			cp3 = 1;
		  }
		  }


		  int isedge1 = 0;
		  int isedge2 = 0;
		  int isedge3 = 0;
		  if (dim == 3)
		  {
		  INDEX_2 i2 = new INDEX_2();
		  i2 = INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
		  isedge1 = edges.Used(i2);
		  i2.Sort();
		  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
		  {
			  isedge1 = 1;
			  ep1 = 1;
			  ep2 = 1;
		  }

		  i2 = INDEX_2(el.PNumMod(j + 1), el.PNumMod(j + 2));
		  isedge2 = edges.Used(i2);
		  i2.Sort();
		  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
		  {
			  isedge2 = 1;
			  ep2 = 1;
			  ep3 = 1;
		  }
		  i2 = INDEX_2(el.PNumMod(j + 2), el.PNumMod(j + 3));
		  isedge3 = edges.Used(i2);
		  i2.Sort();
		  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
		  {
			  isedge3 = 1;
			  ep1 = 1;
			  ep3 = 1;
		  }

		  // cout << " isedge " << isedge1 << " \t " << isedge2 << " \t " << isedge3 << endl;  

		  if (!sing_face)
		  {
				  /*
				    if (!isedge1)  { cp1 |= ep1; cp2 |= ep2; }
				    if (!isedge2)  { cp2 |= ep2; cp3 |= ep3; }
				    if (!isedge3)  { cp3 |= ep3; cp1 |= ep1; }
				  */
				  ep1 |= facepoint [el.PNumMod(j)] != 0;
				  ep2 |= facepoint [el.PNumMod(j + 1)] != 0;
				  ep3 |= facepoint [el.PNumMod(j + 2)] != 0;


				  isedge1 |= face_edges.Used(INDEX_2.Sort(el.PNumMod(j), el.PNumMod(j + 1)));
				  isedge2 |= face_edges.Used(INDEX_2.Sort(el.PNumMod(j + 1), el.PNumMod(j + 2)));
				  isedge3 |= face_edges.Used(INDEX_2.Sort(el.PNumMod(j + 2), el.PNumMod(j + 3)));
		  }
		  }

		  if (dim == 2)
		  {
		  INDEX_2 i2 = new INDEX_2();
		  i2 = INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
		  i2.Sort();
		  isedge1 = edges.Used(i2);
		  if (isedge1 != 0)
		  {
			  ep1 = 1;
			  ep2 = 1;
		  }

		  i2 = INDEX_2(el.PNumMod(j + 1), el.PNumMod(j + 2));
		  i2.Sort();
		  isedge2 = edges.Used(i2);
		  if (isedge2 != 0)
		  {
			  ep2 = 1;
			  ep3 = 1;
		  }
		  i2 = INDEX_2(el.PNumMod(j + 2), el.PNumMod(j + 3));
		  i2.Sort();
		  isedge3 = edges.Used(i2);
		  if (isedge3 != 0)
		  {
			  ep1 = 1;
			  ep3 = 1;
		  }


		  }


		  /*
		    cout << " used " << face_edges.Used (INDEX_2::Sort (el.PNumMod(j), el.PNumMod(j+1))) << endl; 
	
		    cout << " isedge " << isedge1 << " \t " << isedge2 << " \t " << isedge3 << endl; 
		    cout << " ep " << ep1 << "\t" << ep2 << " \t " << ep3 << endl; 
		    cout << " cp " << cp1 << "\t" << cp2 << " \t " << cp3 << endl; 
		  */



		  if (isedge1 + isedge2 + isedge3 == 0)
		  {
		  if (ep1 == 0 && ep2 == 0 && ep3 == 0)
		  {
			type = HP_TRIG;
		  }

		  if (ep1 != 0 && ep2 == 0 && ep3 == 0)
		  {
			type = HP_TRIG_SINGCORNER;
		  }

		  if (ep1 != 0 && ep2 != 0 && ep3 == 0)
		  {
			type = HP_TRIG_SINGCORNER12;
		  }

		  if (ep1 != 0 && ep2 != 0 && ep3 != 0)
		  {
			  if (dim == 2)
			  {
					type = HP_TRIG_SINGCORNER123_2D;
			  }
			  else
			  {
			type = HP_TRIG_SINGCORNER123;
			  }
		  }

		  if (type != HP_NONE)
		  {
			  pnums[0] = el.PNumMod(j);
			  pnums[1] = el.PNumMod(j + 1);
			  pnums[2] = el.PNumMod(j + 2);
			  break;
		  }
		  }

		  if (isedge1 != 0 && isedge2 == 0 && isedge3 == 0)
		  {
		  int code = 0;
		  if (cp1 != 0)
		  {
			  code += 1;
		  }
		  if (cp2 != 0)
		  {
			  code += 2;
		  }
		  if (ep3 != 0)
		  {
			  code += 4;
		  }

		  HPREF_ELEMENT_TYPE[] types = {HP_TRIG_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER12, HP_TRIG_SINGEDGECORNER3, HP_TRIG_SINGEDGECORNER13, HP_TRIG_SINGEDGECORNER23, HP_TRIG_SINGEDGECORNER123};
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: type = types[code];
		  type.CopyFrom(types[code]);
		  pnums[0] = el.PNumMod(j);
		  pnums[1] = el.PNumMod(j + 1);
		  pnums[2] = el.PNumMod(j + 2);
		  break;
		  }


		  if (isedge1 != 0 && isedge2 == 0 && isedge3 != 0)
		  {
		  if (cp3 == 0)
		  {
			  if (cp2 == 0)
			  {
				  type = HP_TRIG_SINGEDGES;
			  }
			  else
			  {
				  type = HP_TRIG_SINGEDGES2;
			  }
		  }
		  else
		  {
			  if (cp2 == 0)
			  {
				  type = HP_TRIG_SINGEDGES3;
			  }
			  else
			  {
				  type = HP_TRIG_SINGEDGES23;
			  }
		  }

		  pnums[0] = el.PNumMod(j);
		  pnums[1] = el.PNumMod(j + 1);
		  pnums[2] = el.PNumMod(j + 2);
		  break;
		  }

		  if (isedge1 != 0 && isedge2 != 0 && isedge3 != 0)
		  {
		  type = HP_TRIG_3SINGEDGES;
		  pnums[0] = el.PNumMod(j);
		  pnums[1] = el.PNumMod(j + 1);
		  pnums[2] = el.PNumMod(j + 2);
		  break;
		  }
	  }

	  for (int k = 0;k < 3;k++)
	  {
		  el[k] = pnums[k];
	  }
	  /*if(type != HP_NONE) 
	    {
	     
	    cout << " TRIG with pnums " << pnums[0] << "\t"  << 
	    pnums[1] << "\t"  << pnums[2] << endl; 
	    cout << " type "  << type << endl; 
	    }
	  */
	  return (type);
	}
	#endif
	public static HPREF_ELEMENT_TYPE ClassifyQuad(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint, int dim, FaceDescriptor fd)
	{
	  HPREF_ELEMENT_TYPE type = HP_NONE;

	  int ep1 = -1;
	  int ep2 = -1;
	  int ep3 = -1;
	  int ep4 = -1;
	  int cp1 = -1;
	  int cp2 = -1;
	  int cp3 = -1;
	  int cp4 = -1;
	  int isedge1;
	  int isedge2;
	  int isedge3;
	  int isedge4;

	  *testout << "edges = " << edges << "\n";

	  for (int j = 1; j <= 4; j++)
	  {
		  ep1 = edgepoint.Test(el.PNumMod(j));
		  ep2 = edgepoint.Test(el.PNumMod(j + 1));
		  ep3 = edgepoint.Test(el.PNumMod(j + 2));
		  ep4 = edgepoint.Test(el.PNumMod(j + 3));

		  if (dim == 2)
		  {
			  ep1 = edgepoint_dom.Used(INDEX_2(el.GetIndex(), el.PNumMod(j)));
			  ep2 = edgepoint_dom.Used(INDEX_2(el.GetIndex(), el.PNumMod(j + 1)));
			  ep3 = edgepoint_dom.Used(INDEX_2(el.GetIndex(), el.PNumMod(j + 2)));
			  ep4 = edgepoint_dom.Used(INDEX_2(el.GetIndex(), el.PNumMod(j + 3)));
		  }

		  cp1 = cornerpoint.Test(el.PNumMod(j));
		  cp2 = cornerpoint.Test(el.PNumMod(j + 1));
		  cp3 = cornerpoint.Test(el.PNumMod(j + 2));
		  cp4 = cornerpoint.Test(el.PNumMod(j + 3));

		  ep1 |= cp1;
		  ep2 |= cp2;
		  ep3 |= cp3;
		  ep4 |= cp4;

		  int[] p = {el.PNumMod(j), el.PNumMod(j + 1), el.PNumMod(j + 2), el.PNumMod(j + 4)};
		  //int epp[4] = { ep1, ep2, ep3, ep4}; 
		  int[] cpp = {cp1, cp2, cp3, cp4};
		  for (int k = 0;k < 0;k++)
		  {
			  INDEX_2 i2a = INDEX_2.Sort(p[k], p[(k + 1) % 4]);
			  INDEX_2 i2b = INDEX_2.Sort(p[k], p[(k - 1) % 4]);
			  if (!edges.Used(i2a) && !edges.Used(i2b))
			  {
				cpp[k] = 1;
			  }
		  }
		  cp1 = cpp[0];
		  cp2 = cpp[1];
		  cp3 = cpp[2];
		  cp4 = cpp[3];


		  if (dim == 3)
		  {
			  INDEX_2 i2 = new INDEX_2();
			  i2 = INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
			  // i2.Sort();
			  isedge1 = edges.Used(i2);
			  i2.Sort();
			  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
			  {
				  isedge1 = 1;
				  ep1 = 1;
				  ep2 = 1;
			  }
			  i2 = INDEX_2(el.PNumMod(j + 1), el.PNumMod(j + 2));
			  // i2.Sort();
			  isedge2 = edges.Used(i2);
			  i2.Sort();
			  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
			  {
				  isedge2 = 1;
				  ep2 = 1;
				  ep3 = 1;
			  }
			  i2 = INDEX_2(el.PNumMod(j + 2), el.PNumMod(j + 3));
			  // i2.Sort();
			  isedge3 = edges.Used(i2);
			  i2.Sort();
			  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
			  {
				  isedge3 = 1;
				  ep3 = 1;
				  ep4 = 1;
			  }
			  i2 = INDEX_2(el.PNumMod(j + 3), el.PNumMod(j + 4));
			  // i2.Sort();
			  isedge4 = edges.Used(i2);
			  i2.Sort();
			  if (surf_edges.Used(i2) && surf_edges.Get(i2) != fd.SurfNr() + 1 && (face_edges.Get(i2) == -1 || face_edges.Get(i2) == fd.DomainIn() || face_edges.Get(i2) == fd.DomainOut()))
			  {
				  isedge4 = 1;
				  ep4 = 1;
				  ep1 = 1;
			  }


			  //MH***********************************************************************************************************
			  if (ep1 != 0)
			  {
				if (edgepoint.Test(p[0]))
				{
					INDEX_2 i2a = INDEX_2.Sort(p[0], p[1]);
					INDEX_2 i2b = INDEX_2.Sort(p[0], p[3]);
					if (!edges.Used(i2a) && !edges.Used(i2b))
					{
					  cp1 = 1;
					}
				}
			  }
			  if (ep2 != 0)
			  {
				if (edgepoint.Test(p[1]))
				{
					INDEX_2 i2a = INDEX_2.Sort(p[0], p[1]);
					INDEX_2 i2b = INDEX_2.Sort(p[1], p[2]);
					if (!edges.Used(i2a) && !edges.Used(i2b))
					{
					  cp2 = 1;
					}
				}
			  }
			  if (ep3 != 0)
			  {
				if (edgepoint.Test(p[2]))
				{
					INDEX_2 i2a = INDEX_2.Sort(p[2], p[1]);
					INDEX_2 i2b = INDEX_2.Sort(p[3], p[2]);
					if (!edges.Used(i2a) && !edges.Used(i2b))
					{
					  cp3 = 1;
					}
				}
			  }
			  if (ep4 != 0)
			  {
				if (edgepoint.Test(p[3]))
				{
					INDEX_2 i2a = INDEX_2.Sort(p[0], p[3]);
					INDEX_2 i2b = INDEX_2.Sort(p[3], p[2]);
					if (!edges.Used(i2a) && !edges.Used(i2b))
					{
					  cp4 = 1;
					}
				}
			  }
			  //MH*****************************************************************************************************************************
		  }
		  else
		  {
			  INDEX_2 i2 = new INDEX_2();
			  i2 = INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
			  i2.Sort();
			  isedge1 = edges.Used(i2);
			  if (isedge1 != 0)
			  {
				  ep1 = 1;
				  ep2 = 1;
			  }
			  i2 = INDEX_2(el.PNumMod(j + 1), el.PNumMod(j + 2));
			  i2.Sort();
			  isedge2 = edges.Used(i2);
			  if (isedge2 != 0)
			  {
				  ep2 = 1;
				  ep3 = 1;
			  }
			  i2 = INDEX_2(el.PNumMod(j + 2), el.PNumMod(j + 3));
			  i2.Sort();
			  isedge3 = edges.Used(i2);

			  if (isedge3 != 0)
			  {
				  ep3 = 1;
				  ep4 = 1;
			  }
			  i2 = INDEX_2(el.PNumMod(j + 3), el.PNumMod(j + 4));
			  i2.Sort();
			  isedge4 = edges.Used(i2);
			  if (isedge4 != 0)
			  {
				  ep4 = 1;
				  ep1 = 1;
			  }
		  }

		  int sumcp = cp1 + cp2 + cp3 + cp4;
		  int sumep = ep1 + ep2 + ep3 + ep4;
		  int sumedge = isedge1 + isedge2 + isedge3 + isedge4;

		  *testout << "isedge = " << isedge1 << isedge2 << isedge3 << isedge4 << "\n";
		  *testout << "iscp = " << cp1 << cp2 << cp3 << cp4 << "\n";
		  *testout << "isep = " << ep1 << ep2 << ep3 << ep4 << "\n";

		  switch (sumedge)
		  {
			case 0:
			{
				switch (sumep)
				{
				  case 0:
					type = HP_QUAD;
					break;
				  case 1:
					if (ep1 != 0)
					{
						type = HP_QUAD_SINGCORNER;
					}
					break;
				  case 2:
				  {
					  if (ep1 != 0 && ep2 != 0)
					  {
						  type = HP_QUAD_0E_2VA;
					  }
					  if (ep1 != 0 && ep3 != 0)
					  {
						  type = HP_QUAD_0E_2VB;
					  }
					  break;
				  }
				  case 3:
					if (ep4 == 0)
					{
						type = HP_QUAD_0E_3V;
					}
					break;
				  case 4:
					type = HP_QUAD_0E_4V;
					break;
				}
				break;
			}
			case 1:
			{
				if (isedge1 != 0)
				{
					switch (cp1 + cp2 + ep3 + ep4)
					{
					  case 0:
						type = HP_QUAD_SINGEDGE;
						break;
					  case 1:
					  {
						  if (cp1 != 0)
						  {
							  type = HP_QUAD_1E_1VA;
						  }
						  if (cp2 != 0)
						  {
							  type = HP_QUAD_1E_1VB;
						  }
						  if (ep3 != 0)
						  {
							  type = HP_QUAD_1E_1VC;
						  }
						  if (ep4 != 0)
						  {
							  type = HP_QUAD_1E_1VD;
						  }
						  break;
					  }
					  case 2:
					  {
						  if (cp1 != 0 && cp2 != 0)
						  {
							  type = HP_QUAD_1E_2VA;
						  }
						  if (cp1 != 0 && ep3 != 0)
						  {
							  type = HP_QUAD_1E_2VB;
						  }
						  if (cp1 != 0 && ep4 != 0)
						  {
							  type = HP_QUAD_1E_2VC;
						  }
						  if (cp2 != 0 && ep3 != 0)
						  {
							  type = HP_QUAD_1E_2VD;
						  }
						  if (cp2 != 0 && ep4 != 0)
						  {
							  type = HP_QUAD_1E_2VE;
						  }
						  if (ep3 != 0 && ep4 != 0)
						  {
							  type = HP_QUAD_1E_2VF;
						  }
						  break;
					  }
					  case 3:
					  {
						  if (cp1 != 0 && cp2 != 0 && ep3 != 0)
						  {
							  type = HP_QUAD_1E_3VA;
						  }
						  if (cp1 != 0 && cp2 != 0 && ep4 != 0)
						  {
							  type = HP_QUAD_1E_3VB;
						  }
						  if (cp1 != 0 && ep3 != 0 && ep4 != 0)
						  {
							  type = HP_QUAD_1E_3VC;
						  }
						  if (cp2 != 0 && ep3 != 0 && ep4 != 0)
						  {
							  type = HP_QUAD_1E_3VD;
						  }
						  break;
					  }
					  case 4:
					  {
						  type = HP_QUAD_1E_4V;
						  break;
					  }
					}
				}
				break;
			}
			case 2:
			{
				if (isedge1 != 0 && isedge4 != 0)
				{
					if (cp2 == 0 && ep3 == 0 && cp4 == 0)
					{
					  type = HP_QUAD_2E;
					}

					if (cp2 != 0 && ep3 == 0 && cp4 == 0)
					{
					  type = HP_QUAD_2E_1VA;
					}
					if (cp2 == 0 && ep3 != 0 && cp4 == 0)
					{
					  type = HP_QUAD_2E_1VB;
					}
					if (cp2 == 0 && ep3 == 0 && cp4 != 0)
					{
					  type = HP_QUAD_2E_1VC;
					}

					if (cp2 != 0 && ep3 != 0 && cp4 == 0)
					{
					  type = HP_QUAD_2E_2VA;
					}
					if (cp2 != 0 && ep3 == 0 && cp4 != 0)
					{
					  type = HP_QUAD_2E_2VB;
					}
					if (cp2 == 0 && ep3 != 0 && cp4 != 0)
					{
					  type = HP_QUAD_2E_2VC;
					}

					if (cp2 != 0 && ep3 != 0 && cp4 != 0)
					{
					  type = HP_QUAD_2E_3V;
					}
				}
				if (isedge1 != 0 && isedge3 != 0)
				{
					switch (sumcp)
					{
					  case 0:
						type = HP_QUAD_2EB_0V;
						break;
					  case 1:
					  {
						  if (cp1 != 0)
						  {
							  type = HP_QUAD_2EB_1VA;
						  }
						  if (cp2 != 0)
						  {
							  type = HP_QUAD_2EB_1VB;
						  }
						  break;
					  }
					  case 2:
					  {
						  if (cp1 != 0 && cp2 != 0)
						  {
							  type = HP_QUAD_2EB_2VA;
						  }
						  if (cp1 != 0 && cp3 != 0)
						  {
							  type = HP_QUAD_2EB_2VB;
						  }
						  if (cp1 != 0 && cp4 != 0)
						  {
							  type = HP_QUAD_2EB_2VC;
						  }
						  if (cp2 != 0 && cp4 != 0)
						  {
							  type = HP_QUAD_2EB_2VD;
						  }
						  break;
					  }
					  case 3:
					  {
						  if (cp1 != 0 && cp2 != 0 && cp3 != 0)
						  {
							  type = HP_QUAD_2EB_3VA;
						  }
						  if (cp1 != 0 && cp2 != 0 && cp4 != 0)
						  {
							  type = HP_QUAD_2EB_3VB;
						  }
						  break;
					  }
					  case 4:
					  {
						  type = HP_QUAD_2EB_4V;
						  break;
					  }
					}
				}
				break;
			}

			case 3:
			{
				if (isedge1 != 0 && isedge2 != 0 && isedge4 != 0)
				{
					if (cp3 == 0 && cp4 == 0)
					{
						type = HP_QUAD_3E;
					}
					if (cp3 != 0 && cp4 == 0)
					{
						type = HP_QUAD_3E_3VA;
					}
					if (cp3 == 0 && cp4 != 0)
					{
						type = HP_QUAD_3E_3VB;
					}
					if (cp3 != 0 && cp4 != 0)
					{
						type = HP_QUAD_3E_4V;
					}
				}
				break;
			}

			case 4:
			{
				type = HP_QUAD_4E;
				break;
			}
		  }

		  if (type != HP_NONE)
		  {
			  int[] pnums = new int[4];
			  pnums[0] = el.PNumMod(j);
			  pnums[1] = el.PNumMod(j + 1);
			  pnums[2] = el.PNumMod(j + 2);
			  pnums[3] = el.PNumMod(j + 3);
			  for (int k = 0;k < 4;k++)
			  {
				  el[k] = pnums[k];
			  }

			  /*  cout << " QUAD with pnums " << pnums[0] << "\t"  << 
			      pnums[1] << "\t"  << pnums[2] << "\t"  << pnums[3] 
			      << endl << " of type " << type << endl; */

			  break;
		  }
	  }
	  if (type == HP_NONE)
	  {
		  (*testout) << "undefined element" << "\n" << "cp = " << cp1 << cp2 << cp3 << cp4 << "\n" << "ep = " << ep1 << ep2 << ep3 << ep4 << "\n" << "isedge = " << isedge1 << isedge2 << isedge3 << isedge4 << "\n";
	  }

	  *testout << "quad type = " << type << "\n";

	  return new HPREF_ELEMENT_TYPE(type);
	}


	public static HPREF_ELEMENT_TYPE ClassifyHex(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint)
	{
	  HPREF_ELEMENT_TYPE type = HP_NONE;

	  // implementation only for HP_HEX_1F_0E_0V
	  //                         HP_HEX_1FA_1FB_0E_0V
	  //                         HP_HEX 
	  // up to now other cases are refine dummies 

	  // indices of bot,top-faces combinations
	  int[][] index =
	  {
		  new int[] {0, 1},
		  new int[] {1, 0},
		  new int[] {2, 4},
		  new int[] {4, 2},
		  new int[] {3, 5},
		  new int[] {5, 3}
	  };
	  int[] p = new int[8];
	  ELEMENT_FACE[] elfaces = MeshTopology.GetFaces1(HEX);
	  ELEMENT_EDGE[] eledges = MeshTopology.GetEdges1(HEX);

	  for (int m = 0;m < 6 && type == HP_NONE;m++)
	  {
		for (int j = 0;j < 4 && type == HP_NONE;j++)
		{
		int[] point_sing = {0, 0, 0, 0, 0, 0, 0, 0};
		int[] face_sing = {0, 0, 0, 0, 0, 0};
		int[] edge_sing = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		int spoint = 0;
		int sface = 0;
		int sedge = 0;
		for (int l = 0;l < 4;l++)
		{
			p[l] = elfaces[index[m][0]][(4 - j - l) % 4];
			p[l + 4] = elfaces[index[m][1]][(j + l) % 4];
		}

		for (int l = 0;l < 8;l++)
		{
		  if (cornerpoint.Test(el.PNum(p[l])))
		  {
			  point_sing[p[l] - 1] = 3;
			  spoint++;
		  }
		  else if (edgepoint.Test(el.PNum(p[l])))
		  {
			  point_sing[p[l] - 1] = 2;
		  }
		  else if (facepoint[el.PNum(p[l])] == -1 || facepoint[el.PNum(p[l])] == el.GetIndex())
		  {
			point_sing[p[l] - 1] = 1;
		  }
		}

		for (int k = 0;k < 12;k++)
		{
				INDEX_2 i2 = INDEX_2.Sort(el.PNum(p[eledges[k][0] - 1]),el.PNum(p[eledges[k][1] - 1]));
				if (edges.Used(i2))
				{
					edge_sing[k] = 2;
					sedge++;
				}
				else
				{
					edge_sing[k] = face_edges.Used(i2);
				}
		}

		for (int k = 0;k < 6;k++)
		{
				INDEX_3 i3 = new INDEX_3();


				INDEX_4 i4 = new INDEX_4(el.pnums[p[elfaces[k][0] - 1] - 1], el.pnums[p[elfaces[k][1] - 1] - 1], el.pnums[p[elfaces[k][2] - 1] - 1],el.pnums[p[elfaces[k][3] - 1] - 1]);
				i4.Sort();
				i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3());

				if (faces.Used(i3))
				{

					int domnr = faces.Get(i3);
					if (domnr == -1 || domnr == el.GetIndex())
					{
						face_sing[k] = 1;
						sface++;
					}

				}
		}

		if (sface == 0 && sedge == 0 && spoint == 0)
		{
			type = HP_HEX;
		}
		if (sedge == 0 && spoint == 0)
		{
			if (face_sing[0] != 0 && face_sing[2] != 0 && sface == 2)
			{
			  type = HP_HEX_1FA_1FB_0E_0V;
			}
			if (face_sing[0] != 0 && sface == 1)
			{
			  type = HP_HEX_1F_0E_0V;
			}
		}

		el.type = type;

		if (type != HP_NONE)
		{
			int[] pnums = new int[8];
			for (int l = 0;l < 8;l++)
			{
				pnums[l] = el[p[l] - 1];
			}
			for (int l = 0;l < 8;l++)
			{
				el[l] = pnums[l];
			}
			/* cout << " HEX with pnums " << pnums[0] << "\t"  << 
				   pnums[1] << "\t"  << pnums[2] << "\t"  << pnums[3] << "\t"  << 
				   pnums[4] << "\t"  <<  pnums[5] << endl << " of type " << type << endl; */
			break;
		}
		}
	  }

	  return (type);

	}

	public static HPREF_ELEMENT_TYPE ClassifySegm(HPRefElement hpel, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint)
	{

	  int cp1 = cornerpoint.Test(hpel[0]);
	  int cp2 = cornerpoint.Test(hpel[1]);

	  INDEX_2 i2 = new INDEX_2();
	  i2 = INDEX_2(hpel[0], hpel[1]);
	  i2.Sort();
	  if (!edges.Used(i2))
	  {
		  cp1 = edgepoint.Test(hpel[0]);
		  cp2 = edgepoint.Test(hpel[1]);
	  }

	  if (!edges.Used(i2) && !face_edges.Used(i2))
	  {
		  if (facepoint[hpel[0]] != 0)
		  {
			  cp1 = 1;
		  }
		  if (facepoint[hpel[1]] != 0)
		  {
			  cp2 = 1;
		  }
	  }

	  if (edges.Used(i2) && !face_edges.Used(i2))
	  {
		  if (facepoint[hpel[0]])
		  {
			  cp1 = 1;
		  }
		  if (facepoint[hpel[1]])
		  {
			  cp2 = 1;
		  }
	  }

	  if (cp1 == 0 && cp2 == 0)
	  {
		hpel.type = HP_SEGM;
	  }
	  else if (cp1 && !cp2)
	  {
		hpel.type = HP_SEGM_SINGCORNERL;
	  }
	  else if (!cp1 && cp2)
	  {
		hpel.type = HP_SEGM_SINGCORNERR;
	  }
	  else
	  {
		hpel.type = HP_SEGM_SINGCORNERS;
	  }

	  // cout << " SEGM found with " << hpel[0] << " \t" << hpel[1] << endl << " of type " << hpel.type << endl; 
	  return (hpel.type);
	}


	public static HPREF_ELEMENT_TYPE ClassifyPyramid(HPRefElement el, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, Array<int, PointIndex.BASE> facepoint)
	{
	  HPREF_ELEMENT_TYPE type = HP_NONE;

	  // implementation only for HP_PYRAMID
	  //                         HP_PYRAMID_0E_1V
	  //                         HP_PYRAMID_EDGES
	  //                         HP_PYRAMID_1FB_0E_1VA
	  // up to now other cases are refine dummies 

	  // indices of bot,top-faces combinations
	  // int index[6][2] = {{0,1},{1,0},{2,4},{4,2},{3,5},{5,3}}; 

	  ELEMENT_FACE[] elfaces = MeshTopology.GetFaces1(PYRAMID);
	  ELEMENT_EDGE[] eledges = MeshTopology.GetEdges1(PYRAMID);

	  int[] point_sing = {0, 0, 0, 0, 0};
	  int[] face_sing = {0, 0, 0, 0, 0};
	  int[] edge_sing = {0, 0, 0, 0, 0, 0, 0, 0};

	  int spoint = 0;
	  int sedge = 0;
	  int sface = 0;

	  for (int m = 0;m < 4 && type == HP_NONE;m++)
	  {
		  int[] p = {m % 4, m % 4 + 1, m % 4 + 2, m % 4 + 3, 4};

		  for (int l = 0;l < 5;l++)
		  {
		  if (cornerpoint.Test(el.pnums[p[l]]))
		  {
			point_sing[l] = 3;
		  }

			  else if (edgepoint.Test(el.pnums[p[l]]))
			  {
			point_sing[l] = 2;
			  }

		  else if (facepoint[el.pnums[p[l]]] == -1 || facepoint[el.pnums[p[l]]] == el.GetIndex())
		  {
			point_sing[l] = 1;
		  }

		  spoint += point_sing[l];
		  }

		  for (int k = 0;k < 8;k++)
		  {
		  INDEX_2 i2 = INDEX_2.Sort(el.pnums[p[eledges[k][0] - 1]], el.pnums[p[eledges[k][1] - 1]]);
		  if (edges.Used(i2))
		  {
			edge_sing[k] = 2;
		  }
		  else
		  {
			edge_sing[k] = face_edges.Used(i2);
		  }

		  sedge += edge_sing[k];
		  }

		  for (int k = 0;k < 5;k++)
		  {
		  INDEX_3 i3 = new INDEX_3();
		  INDEX_4 i4 = new INDEX_4(el.pnums[p[elfaces[k][0] - 1]], el.pnums[p[elfaces[k][1] - 1]], el.pnums[p[elfaces[k][2] - 1]], el.pnums[p[elfaces[k][3] - 1]]);
		  i4.Sort();
		  i3 = INDEX_3(i4.I1(), i4.I2(), i4.I3());

		  if (faces.Used(i3))
		  {

			  int domnr = faces.Get(i3);
			  if (domnr == -1 || domnr == el.GetIndex())
			  {
			face_sing[k] = 1;
			  }
		  }
		  sface += face_sing[k];
		  }

		  if (sface == 0 && spoint == 0 && sedge == 0)
		  {
			  return (HP_PYRAMID);
		  }

		  if (sface == 0 && sedge == 0 && point_sing[p[0]] == spoint)
		  {
		type = HP_PYRAMID_0E_1V;
		  }

		  if (sface == 0 && edge_sing[0] + edge_sing[2] == sedge && spoint == point_sing[0] + point_sing[1] + point_sing[3])
		  {
		type = HP_PYRAMID_EDGES;
		  }

		  if (sface != 0 && sface == face_sing[0] && spoint == point_sing[4] + 2)
		  {
		type = HP_PYRAMID_1FB_0E_1VA;
		  }


		  if (type != HP_NONE)
		  {
		  int[] pnums = new int[8];
		  for (int l = 0;l < 5;l++)
		  {
			  pnums[l] = el[p[l]];
		  }
		  for (int l = 0;l < 5;l++)
		  {
			  el[l] = pnums[l];
		  }
		  el.type = type;
		  break;
		  }
	  }

	  return (type);

	}

	// find inner point



	public static void Minimize(Array<Vec3d> a, Array<double> c, int[] act, ref Vec < 3> x, ref double f, int[] sol)
	{
	  int[] act1 = new int[4];
	  Mat < 3> m, inv;
	  Vec < 3> rs, xmax, center;

	  f = 1e99;

	  for (int j = 0; j < 5; j++)
	  {
		  for (int hk = 0, k = 0; hk < 4; hk++)
		  {
		  if (hk == j)
		  {
			  k++;
		  }
		  act1[hk] = act[k];
		  k++;
		  }

		  for (int k = 0; k < 3; k++)
		  {
		  m(k, 0) = a[act1[0]].X() - a[act1[k + 1]].X();
		  m(k, 1) = a[act1[0]].Y() - a[act1[k + 1]].Y();
		  m(k, 2) = a[act1[0]].Z() - a[act1[k + 1]].Z();
		  rs(k) = c[act1[k + 1]] - c[act1[0]];
		  }

		  /*
		  (*testout) << "act1 = "
			 << act1[0] << " "
			 << act1[1] << " "
			 << act1[2] << " "
			 << act1[3] << endl;
		  (*testout) << "Det = " << Det(m) << endl;
		  */

		  if (Math.Abs(Det(m)) > 1e-10)
		  {
		  CalcInverse(m, inv);
		  xmax = inv * rs;

		  double fmax = -1e10;
		  for (int k = 0; k < 5; k++)
		  {
			  double hd = xmax(0) * a[act[k]].X() + xmax(1) * a[act[k]].Y() + xmax(2) * a[act[k]].Z() + c[act[k]];
			  if (hd > fmax)
			  {
				  fmax = hd;
			  }
		  }

		  if (fmax < f)
		  {
			  f = fmax;
			  x = xmax;
			  for (int k = 0; k < 4; k++)
			  {
			sol[k] = act1[k];
			  }
		  }
		  }
	  }
	}




//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename POINTArray, typename FACEArray>
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: There is no equivalent in C# to templates on variables:
	internal static int FindInnerPoint_timer = NgProfiler.CreateTimer("FindInnerPoint");

	public static int FindInnerPoint(POINTArray points, FACEArray faces, ref Point3d p)
	{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//  static int timer = NgProfiler::CreateTimer("FindInnerPoint");
	  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(FindInnerPoint_timer);

	  Array<Vec3d> a = new Array<Vec3d>();
	  Array<double> c = new Array<double>();
	  Mat < 3> m, inv;
	  Vec < 3> rs, x = 0.0, center;
	  double f;

	  int nf = faces.Size();

	  // minimize_x  max_i  a_i x + c_i

	  a.SetSize(nf + 4);
	  c.SetSize(nf + 4);

	  for (int i = 0; i < nf; i++)
	  {
		  Point3d p1 = points.Get(faces[i][0]);
		  a[i] = Cross(points.Get(faces[i][1]) - p1, points.Get(faces[i][2]) - p1);
		  a[i] /= a[i].Length();
		  c[i] = - (a[i].X() * p1.X() + a[i].Y() * p1.Y() + a[i].Z() * p1.Z());
	  }

	  /*
	  center = 0;
	  for (int i = 0; i < points.Size(); i++)
	    center += Vec<3> (points[i]);
	  center /= points.Size();
	  */

	  center = 0;
	  for (int i = 0; i < faces.Size(); i++)
	  {
		for (int j = 0; j < 3; j++)
		{
		  center += Vec < 3> (points.Get(faces[i][j]));
		}
	  }
	  center /= (3 * faces.Size());


	  // (*testout) << "center = " << center << endl;

	  double hmax = 0;
	  for (int i = 0; i < nf; i++)
	  {
		  // const Element2d & el = faces[i];
		  // (*testout) << "el[" << i << "] = " << el << endl;
		  for (int j = 1; j <= 3; j++)
		  {
		  double hi = Dist(points.Get(faces[i].PNumMod(j)), points.Get(faces[i].PNumMod(j + 1)));
		  if (hi > hmax)
		  {
			  hmax = hi;
		  }
		  }
	  }

	  // (*testout) << "hmax = " << hmax << endl;

	  a[nf] = Vec < 3> (1, 0, 0);
	  c[nf] = -center(0) - hmax;
	  a[nf + 1] = Vec < 3> (0, 1, 0);
	  c[nf + 1] = -center(1) - hmax;
	  a[nf + 2] = Vec < 3> (0, 0, 1);
	  c[nf + 2] = -center(2) - hmax;
	  a[nf + 3] = Vec < 3> (-1, -1, -1);
	  c[nf + 3] = center(0) + center(1) + center(2) - 3 * hmax;

	  /*
	  (*testout) << "findip, a now = " << endl << a << endl;
	  (*testout) << "findip, c now = " << endl << c << endl;
	  */

	  int[] act = {0, nf, nf + 1, nf + 2, nf + 3};
	  int[] sol = new int[4];

	  while (true)
	  {
		  /*
		  (*testout) << "try ";
		  for (int j = 0; j < 5; j++)
		(*testout)  << act[j] << " ";
		  */

		  Minimize(a, c, act, ref x, ref f, sol);

		  /*
		  (*testout) << endl << "sol = ";
		  for (int j = 0; j < 4; j++)
		(*testout)  << sol[j] << " ";
	
		  (*testout) << " fmin = " << f << endl;
		  */
		  for (int j = 0; j < 4; j++)
		  {
			  act[j] = sol[j];
		  }

		  bool found = false;
		  double maxval = f;
		  for (int j = 0; j < nf; j++)
		  {
		  double val = x(0) * a[j].X() + x(1) * a[j].Y() + x(2) * a[j].Z() + c[j];
		  if (val > maxval + hmax * 1e-6)
		  {
			  found = true;
			  maxval = val;
			  act[4] = j;
		  }
		  }

		  // (*testout) << "maxval = " << maxval << endl;
		  if (!found)
		  {
			  break;
		  }
	  }

	  // cout << "converged, f = " << f << endl;

	  p = Point3d(x(0), x(1), x(2));
	  // (*testout) << "findip, f = " << f << ", hmax = " << hmax << endl;
	  return (f < -1e-5 * hmax);
	}




	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: There is no equivalent in C# to templates on variables:
	internal static int FindInnerPoint2_timer = NgProfiler.CreateTimer("FindInnerPoint2");

	public static int FindInnerPoint2(POINTArray points, FACEArray faces, ref Point3d p)
	{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//  static int timer = NgProfiler::CreateTimer("FindInnerPoint2");
	  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(FindInnerPoint2_timer);

	  Array<Vec3d> a = new Array<Vec3d>();
	  Array<double> c = new Array<double>();
	  Mat < 3> m, inv;
	  Vec < 3> rs, x, pmin;

	  int nf = faces.Size();

	  a.SetSize(nf);
	  c.SetSize(nf);

	  for (int i = 0; i < nf; i++)
	  {
		  Point3d p1 = points.Get(faces[i][0]);
		  a[i] = Cross(points.Get(faces[i][1]) - p1, points.Get(faces[i][2]) - p1);
		  a[i] /= a[i].Length();
		  c[i] = - (a[i].X() * p1.X() + a[i].Y() * p1.Y() + a[i].Z() * p1.Z());
	  }


	  x = 0;


	  double hmax = 0;
	  for (int i = 0; i < nf; i++)
	  {
		  Element2d el = faces[i];
		  for (int j = 1; j <= 3; j++)
		  {
		  double hi = Dist(points.Get(el.PNumMod(j)), points.Get(el.PNumMod(j + 1)));
		  if (hi > hmax)
		  {
			  hmax = hi;
		  }
		  }
	  }

	  double fmin = 0;

	  for (int i1 = 1; i1 <= nf; i1++)
	  {
		for (int i2 = i1 + 1; i2 <= nf; i2++)
		{
		  for (int i3 = i2 + 1; i3 <= nf; i3++)
		  {
			for (int i4 = i3 + 1; i4 <= nf; i4++)
			{
			m(0, 0) = a.Get(i1).X() - a.Get(i2).X();
			m(0, 1) = a.Get(i1).Y() - a.Get(i2).Y();
			m(0, 2) = a.Get(i1).Z() - a.Get(i2).Z();
			rs(0) = c.Get(i2) - c.Get(i1);

			m(1, 0) = a.Get(i1).X() - a.Get(i3).X();
			m(1, 1) = a.Get(i1).Y() - a.Get(i3).Y();
			m(1, 2) = a.Get(i1).Z() - a.Get(i3).Z();
			rs(1) = c.Get(i3) - c.Get(i1);

			m(2, 0) = a.Get(i1).X() - a.Get(i4).X();
			m(2, 1) = a.Get(i1).Y() - a.Get(i4).Y();
			m(2, 2) = a.Get(i1).Z() - a.Get(i4).Z();
			rs(2) = c.Get(i4) - c.Get(i1);


			if (Math.Abs(Det(m)) > 1e-10)
			{
			CalcInverse(m, inv);
			x = inv * rs;

			double f = -1e10;
			for (int i = 0; i < nf; i++)
			{
				double hd = x(0) * a[i].X() + x(1) * a[i].Y() + x(2) * a[i].Z() + c[i];
				if (hd > f)
				{
					f = hd;
				}
				if (hd > fmin)
				{
					break;
				}
			}

			if (f < fmin)
			{
				fmin = f;
				pmin = x;
			}
			}
			}
		  }
		}
	  }

	  p = Point3d(pmin(0), pmin(1), pmin(2));
	  (*testout) << "fmin = " << fmin << "\n";
	  return (fmin < -1e-3 * hmax);
	}


	// SZ 

	// HP_HEX  ... no refinement
	public static int[][] refhex_splitedges =
	{
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refhex_newelstypes = {HP_HEX, HP_NONE};
	public static int[][] refhex_newels =
	{
		new int[] {1, 2, 3, 4, 5, 6, 7, 8}
	};
	public static HPRef_Struct refhex = new HPRef_Struct(HP_HEX, refhex_splitedges, 0, 0, refhex_newelstypes, refhex_newels);

	// HP_HEX_1F  ... face (1 - 4 - 3 -2) singular 
	public static int[][] refhex_1f_0e_0v_splitedges =
	{
		new int[] {1, 5, 9},
		new int[] {2, 6, 10},
		new int[] {3, 7, 11},
		new int[] {4, 8, 12},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refhex_1f_0e_0v_newelstypes = {HP_HEX, HP_HEX_1F_0E_0V, HP_NONE};
	public static int[][] refhex_1f_0e_0v_newels =
	{
		new int[] {9, 10, 11, 12, 5, 6, 7, 8},
		new int[] {1, 2, 3, 4, 9, 10, 11, 12}
	};
	public static HPRef_Struct refhex_1f_0e_0v = new HPRef_Struct(HP_HEX, refhex_1f_0e_0v_splitedges, 0, 0, refhex_1f_0e_0v_newelstypes, refhex_1f_0e_0v_newels);



	// HP_HEX_1FA_1FB  ... face (1 - 4 - 3 -2) and face (1-2-6-5) singular 
	public static int[][] refhex_1fa_1fb_0e_0v_splitedges =
	{
		new int[] {1, 5, 9},
		new int[] {2, 6, 10},
		new int[] {3, 7, 11},
		new int[] {4, 8, 12},
		new int[] {1, 4, 13},
		new int[] {2, 3, 14},
		new int[] {6, 7, 15},
		new int[] {5, 8, 16},
		new int[] {0, 0, 0}
	};

	public static int[][] refhex_1fa_1fb_0e_0v_splitfaces =
	{
		new int[] {2, 3, 6, 17},
		new int[] {1, 4, 5, 18},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refhex_1fa_1fb_0e_0v_newelstypes = {HP_HEX, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_NONE};
	public static int[][] refhex_1fa_1fb_0e_0v_newels =
	{
		new int[] {18, 17, 11, 12, 16, 15, 7, 8},
		new int[] {13, 14, 3, 4, 18, 17, 11, 12},
		new int[] {5, 6, 10, 9, 16, 15, 17, 18},
		new int[] {1, 2, 14, 13, 9, 10, 17, 18}
	};
	public static HPRef_Struct refhex_1fa_1fb_0e_0v = new HPRef_Struct(HP_HEX, refhex_1fa_1fb_0e_0v_splitedges, refhex_1fa_1fb_0e_0v_splitfaces, 0, refhex_1fa_1fb_0e_0v_newelstypes, refhex_1fa_1fb_0e_0v_newels);



	// Refine Dummies 
	  // HP_HEX_0E_1V
	  public static int[][] refhex_0e_1v_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refhex_0e_1v_newelstypes = {HP_TET_0E_1V, HP_TET, HP_TET, HP_TET, HP_TET, HP_TET, HP_NONE};
	  public static int[][] refhex_0e_1v_newels =
	  {
		  new int[] {1, 5, 2, 4, 0, 0, 0, 0},
		  new int[] {7, 3, 6, 8, 0, 0, 0, 0},
		  new int[] {2, 8, 5, 6, 0, 0, 0, 0},
		  new int[] {2, 8, 6, 3, 0, 0, 0, 0},
		  new int[] {2, 8, 3, 4, 0, 0, 0, 0},
		  new int[] {2, 8, 4, 5, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refhex_0e_1v = new HPRef_Struct(HP_HEX, refhex_0e_1v_splitedges, 0, 0, refhex_0e_1v_newelstypes, refhex_0e_1v_newels);



	// Refine Dummies 
	  // HP_HEX_1E_1V
	  public static int[][] refhex_1e_1v_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refhex_1e_1v_newelstypes = {HP_TET_1E_1VA, HP_TET, HP_TET_0E_1V, HP_TET_0E_1V, HP_TET_0E_1V, HP_TET_0E_1V, HP_NONE};
	  public static int[][] refhex_1e_1v_newels =
	  {
		  new int[] {1, 2, 4, 5, 0, 0, 0, 0},
		  new int[] {7, 3, 6, 8, 0, 0, 0, 0},
		  new int[] {2, 8, 5, 6, 0, 0, 0, 0},
		  new int[] {2, 8, 6, 3, 0, 0, 0, 0},
		  new int[] {2, 8, 3, 4, 0, 0, 0, 0},
		  new int[] {2, 8, 4, 5, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refhex_1e_1v = new HPRef_Struct(HP_HEX, refhex_1e_1v_splitedges, 0, 0, refhex_1e_1v_newelstypes, refhex_1e_1v_newels);


	// Refine Dummies 
	  // HP_HEX_3E_0V
	  public static int[][] refhex_3e_0v_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refhex_3e_0v_newelstypes = {HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_0E_1V, HP_TET, HP_NONE};
	  public static int[][] refhex_3e_0v_newels =
	  {
		  new int[] {1, 2, 3, 6, 0, 0, 0, 0},
		  new int[] {1, 4, 8, 3, 0, 0, 0, 0},
		  new int[] {1, 5, 6, 8, 0, 0, 0, 0},
		  new int[] {1, 6, 3, 8, 0, 0, 0, 0},
		  new int[] {3, 8, 6, 7, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refhex_3e_0v = new HPRef_Struct(HP_HEX, refhex_3e_0v_splitedges, 0, 0, refhex_3e_0v_newelstypes, refhex_3e_0v_newels);



	// Refine Dummies 
	  // HP_HEX_1E_0V 
	  public static int[][] refhex_1e_0v_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refhex_1e_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM, HP_NONE};
	  public static int[][] refhex_1e_0v_newels =
	  {
		  new int[] {1, 4, 5, 2, 3, 6, 0, 0},
		  new int[] {5, 4, 8, 6, 3, 7, 0, 0}
	  };
	  public static HPRef_Struct refhex_1e_0v = new HPRef_Struct(HP_HEX, refhex_1e_0v_splitedges, 0, 0, refhex_1e_0v_newelstypes, refhex_1e_0v_newels);



	  // HP_PRISM  ... no refinement
	  public static int[][] refprism_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_newelstypes = {HP_PRISM, HP_NONE};
	  public static int[][] refprism_newels =
	  {
		  new int[] {1, 2, 3, 4, 5, 6, 0, 0}
	  };
	  public static HPRef_Struct refprism = new HPRef_Struct(HP_PRISM, refprism_splitedges, 0, 0, refprism_newelstypes, refprism_newels);



	  // HP_PRISM_SINGEDGE  ... vertical edge 1-4 is singular
	  public static int[][] refprism_singedge_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {1, 3, 8},
		  new int[] {4, 5, 9},
		  new int[] {4, 6, 10},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_singedge_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_NONE};
	  public static int[][] refprism_singedge_newels =
	  {
		  new int[] {1, 7, 8, 4, 9, 10, 0, 0},
		  new int[] {3, 8, 7, 2, 6, 10, 9, 5}
	  };
	  public static HPRef_Struct refprism_singedge = new HPRef_Struct(HP_PRISM, refprism_singedge_splitedges, 0, 0, refprism_singedge_newelstypes, refprism_singedge_newels);






	  // HP_PRISM_SINGEDGE_V12  vertical edges 1-4 and 2-5 are singular 
	  public static int[][] refprism_singedge_v12_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {1, 3, 8},
		  new int[] {2, 1, 9},
		  new int[] {2, 3, 10},
		  new int[] {4, 5, 11},
		  new int[] {4, 6, 12},
		  new int[] {5, 4, 13},
		  new int[] {5, 6, 14},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_singedge_v12_newelstypes = {HP_HEX, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM, HP_NONE};
	  public static int[][] refprism_singedge_v12_newels =
	  {
		  new int[] {7, 9, 10, 8, 11, 13, 14, 12},
		  new int[] {1, 7, 8, 4, 11, 12, 0, 0},
		  new int[] {2, 10, 9, 5, 14, 13, 0, 0},
		  new int[] {3, 8, 10, 6, 12, 14, 0, 0}
	  };
	  public static HPRef_Struct refprism_singedge_v12 = new HPRef_Struct(HP_PRISM, refprism_singedge_v12_splitedges, 0, 0, refprism_singedge_v12_newelstypes, refprism_singedge_v12_newels);






	  // HP_PRISM_SINGEDGE_H12
	  public static int[][] refprism_singedge_h12_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 1, 8},
		  new int[] {2, 3, 9},
		  new int[] {3, 1, 10},
		  new int[] {4, 6, 12},
		  new int[] {5, 4, 13},
		  new int[] {5, 6, 14},
		  new int[] {6, 4, 15},
		  new int[] {0, 0, 0}
	  };

	  public static int[][] refprism_singedge_h12_splitfaces =
	  {
		  new int[] {2, 1, 3, 11},
		  new int[] {5, 4, 6, 16},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_singedge_h12_newelstypes = {HP_HEX, HP_HEX, HP_PRISM, HP_PRISM, HP_PRISM, HP_NONE};
	  public static int[][] refprism_singedge_h12_newels =
	  {
		  new int[] {1, 8, 11, 7, 4, 13, 16, 12},
		  new int[] {9, 3, 10, 11, 14, 6, 15, 16},
		  new int[] {7, 11, 10, 12, 16, 15, 0, 0},
		  new int[] {2, 9, 11, 5, 14, 16, 0, 0},
		  new int[] {8, 2, 11, 13, 5, 16, 0, 0}
	  };
	  public static HPRef_Struct refprism_singedge_h12 = new HPRef_Struct(HP_PRISM, refprism_singedge_h12_splitedges, refprism_singedge_h12_splitfaces, 0, refprism_singedge_h12_newelstypes, refprism_singedge_h12_newels);






	  // HP_PRISM_SINGEDGE_H1
	  public static int[][] refprism_singedge_h1_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {4, 6, 9},
		  new int[] {5, 6, 10},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_singedge_h1_newelstypes = {HP_HEX, HP_PRISM, HP_NONE};
	  public static int[][] refprism_singedge_h1_newels =
	  {
		  new int[] {1, 2, 8, 7, 4, 5, 10, 9},
		  new int[] {3, 7, 8, 6, 9, 10, 0, 0}
	  };
	  public static HPRef_Struct refprism_singedge_h1 = new HPRef_Struct(HP_PRISM, refprism_singedge_h1_splitedges, 0, 0, refprism_singedge_h1_newelstypes, refprism_singedge_h1_newels);



	//  HP_PRISM_1FA_0E_0V
	  public static int[][] refprism_1fa_0e_0v_splitedges =
	  {
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_0e_0v_newelstypes = {HP_PRISM, HP_PRISM_1FA_0E_0V, HP_NONE};
	  public static int[][] refprism_1fa_0e_0v_newels =
	  {
		  new int[] {16, 17, 18, 4, 5, 6, 0, 0},
		  new int[] {1, 2, 3, 16, 17, 18, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_0e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_0e_0v_splitedges, 0, 0, refprism_1fa_0e_0v_newelstypes, refprism_1fa_0e_0v_newels);

	//  HP_PRISM_1FA_1E_0V
	  public static int[][] refprism_1fa_1e_0v_splitedges =
	  {
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {1, 2, 7},
		  new int[] {1, 3, 12},
		  new int[] {4, 6, 45},
		  new int[] {4, 5, 40},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_1e_0v_splitfaces =
	  {
		  new int[] {1, 2, 4, 19},
		  new int[] {1, 3, 4, 24},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1e_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_NONE};
	  public static int[][] refprism_1fa_1e_0v_newels =
	  {
		  new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		  new int[] {24, 19, 17, 18, 45, 40, 5, 6},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {7, 2, 3, 12, 19, 17, 18, 24}
	  };
	  public static HPRef_Struct refprism_1fa_1e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1e_0v_splitedges, refprism_1fa_1e_0v_splitfaces, 0, refprism_1fa_1e_0v_newelstypes, refprism_1fa_1e_0v_newels);

	//  HP_PRISM_2FA_1E_0V
	  public static int[][] refprism_2fa_1e_0v_splitedges =
	  {
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {1, 2, 7},
		  new int[] {1, 3, 12},
		  new int[] {4, 6, 45},
		  new int[] {4, 5, 40},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_2fa_1e_0v_splitfaces =
	  {
		  new int[] {1, 2, 4, 19},
		  new int[] {1, 3, 4, 24},
		  new int[] {4, 1, 5, 31},
		  new int[] {4, 1, 6, 36},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1e_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_NONE};
	  public static int[][] refprism_2fa_1e_0v_newels =
	  {
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {24, 19, 17, 18, 36, 31, 29, 30},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {12, 7, 2, 3, 24, 19, 17, 18},
		  new int[] {4, 45, 40, 28, 36, 31, 0, 0},
		  new int[] {40, 45, 6, 5, 31, 36, 30, 29}
	  };
	  public static HPRef_Struct refprism_2fa_1e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1e_0v_splitedges, refprism_2fa_1e_0v_splitfaces, 0, refprism_2fa_1e_0v_newelstypes, refprism_2fa_1e_0v_newels);

	//  HP_PRISM_1FB_0E_0V   ... quad face 1-2-4-5
	  public static int[][] refprism_1fb_0e_0v_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {4, 6, 9},
		  new int[] {5, 6, 10},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fb_0e_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PRISM, HP_NONE};
	  public static int[][] refprism_1fb_0e_0v_newels =
	  {
		  new int[] {1, 4, 5, 2, 7, 9, 10, 8},
		  new int[] {7, 8, 3, 9, 10, 6, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fb_0e_0v = new HPRef_Struct(HP_PRISM, refprism_1fb_0e_0v_splitedges, 0, 0, refprism_1fb_0e_0v_newelstypes, refprism_1fb_0e_0v_newels);


	//  HP_PRISM_1FB_1EA_0V   ... quad face 1-2-4-5
	  public static int[][] refprism_1fb_1ea_0v_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {4, 6, 9},
		  new int[] {5, 6, 10},
		  new int[] {1, 2, 11},
		  new int[] {4, 5, 12},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fb_1ea_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM, HP_NONE};
	  public static int[][] refprism_1fb_1ea_0v_newels =
	  {
		  new int[] {11, 12, 5, 2, 7, 9, 10, 8},
		  new int[] {1, 11, 7, 4, 12, 9, 0, 0},
		  new int[] {7, 8, 3, 9, 10, 6, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fb_1ea_0v = new HPRef_Struct(HP_PRISM, refprism_1fb_1ea_0v_splitedges, 0, 0, refprism_1fb_1ea_0v_newelstypes, refprism_1fb_1ea_0v_newels);

	//  HP_PRISM_1FB_1EC_0V   ... quad face 1-2-4-5 with singular edge 3-6 
	  public static int[][] refprism_1fb_1ec_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fb_1ec_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_HEX_1F_0E_0V, HP_NONE};
	  public static int[][] refprism_1fb_1ec_0v_newels =
	  {
		  new int[] {3, 11, 10, 6, 44, 43, 0, 0},
		  new int[] {12, 9, 10, 11, 45, 42, 43, 44},
		  new int[] {4, 5, 2, 1, 45, 42, 9, 12}
	  };
	  public static HPRef_Struct refprism_1fb_1ec_0v = new HPRef_Struct(HP_PRISM, refprism_1fb_1ec_0v_splitedges, 0, 0, refprism_1fb_1ec_0v_newelstypes, refprism_1fb_1ec_0v_newels);

	//  HP_PRISM_1FA_1FB_1EC_0V   ... bot-trig face, quad face 1-2-4-5 with singular edge 3-6 
	  public static int[][] refprism_1fa_1fb_1ec_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {0, 0, 0}
	  };

	 public static int[][] refprism_1fa_1fb_1ec_0v_splitfaces =
	 {
		 new int[] {2, 3, 5, 21},
		 new int[] {3, 2, 6, 22},
		 new int[] {3, 1, 6, 23},
		 new int[] {1, 3, 4, 24},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_1ec_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_1ec_0v_newels =
	  {
		  new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		  new int[] {24, 21, 22, 23, 45, 42, 43, 44},
		  new int[] {4, 5, 17, 16, 45, 42, 21, 24},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {1, 2, 9, 12, 16, 17, 21, 24}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_1ec_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_1ec_0v_splitedges, refprism_1fa_1fb_1ec_0v_splitfaces, 0, refprism_1fa_1fb_1ec_0v_newelstypes, refprism_1fa_1fb_1ec_0v_newels);


	//  HP_PRISM_1FA_1FB_2EB_0V  
	  public static int[][] refprism_1fa_1fb_2eb_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {4, 5, 40},
		  new int[] {4, 6, 45},
		  new int[] {1, 2, 7},
		  new int[] {0, 0, 0}
	  };

	 public static int[][] refprism_1fa_1fb_2eb_0v_splitfaces =
	 {
		 new int[] {2, 3, 5, 21},
		 new int[] {3, 2, 6, 22},
		 new int[] {3, 1, 6, 23},
		 new int[] {1, 3, 4, 24},
		 new int[] {1, 2, 4, 19},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_2eb_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_2eb_0v_newels =
	  {
		  new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		  new int[] {24, 21, 22, 23, 45, 42, 43, 44},
		  new int[] {40, 5, 17, 19, 45, 42, 21, 24},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {7, 2, 9, 12, 19, 17, 21, 24},
		  new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_2eb_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_2eb_0v_splitedges, refprism_1fa_1fb_2eb_0v_splitfaces, 0, refprism_1fa_1fb_2eb_0v_newelstypes, refprism_1fa_1fb_2eb_0v_newels);

	 //  HP_PRISM_1FA_1FB_2EC_0V 
	  public static int[][] refprism_1fa_1fb_2ec_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {0, 0, 0}
	  };

	 public static int[][] refprism_1fa_1fb_2ec_0v_splitfaces =
	 {
		 new int[] {2, 3, 5, 21},
		 new int[] {3, 2, 6, 22},
		 new int[] {3, 1, 6, 23},
		 new int[] {1, 3, 4, 24},
		 new int[] {2, 1, 5, 20},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_2ec_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_2ec_0v_newels =
	  {
		  new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		  new int[] {24, 21, 22, 23, 45, 42, 43, 44},
		  new int[] {4, 41, 20, 16, 45, 42, 21, 24},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {1, 8, 9, 12, 16, 20, 21, 24},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0},
		  new int[] {5, 41, 42, 17, 20, 21, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_2ec_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_2ec_0v_splitedges, refprism_1fa_1fb_2ec_0v_splitfaces, 0, refprism_1fa_1fb_2ec_0v_newelstypes, refprism_1fa_1fb_2ec_0v_newels);







	//  HP_PRISM_2FA_1FB_1EC_0V   ... trig faces, quad face 1-2-4-5 with singular edge 3-6 
	  public static int[][] refprism_2fa_1fb_1ec_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {0, 0, 0}
	  };

	 public static int[][] refprism_2fa_1fb_1ec_0v_splitfaces =
	 {
		 new int[] {2, 3, 5, 21},
		 new int[] {3, 2, 6, 22},
		 new int[] {3, 1, 6, 23},
		 new int[] {1, 3, 4, 24},
		 new int[] {5, 2, 6, 33},
		 new int[] {6, 5, 3, 34},
		 new int[] {6, 4, 3, 35},
		 new int[] {4, 1, 6, 36},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1fb_1ec_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_NONE};
	  public static int[][] refprism_2fa_1fb_1ec_0v_newels =
	  {
		  new int[] {18, 23, 22, 30, 35, 34, 0, 0},
		  new int[] {24, 21, 22, 23, 36, 33, 34, 35},
		  new int[] {28, 29, 17, 16, 36, 33, 21, 24},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {1, 2, 9, 12, 16, 17, 21, 24},
		  new int[] {6, 43, 44, 30, 34, 35, 0, 0},
		  new int[] {44, 43, 42, 45, 35, 34, 33, 36},
		  new int[] {5, 4, 45, 42, 29, 28, 36, 33}
	  };
	  public static HPRef_Struct refprism_2fa_1fb_1ec_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1fb_1ec_0v_splitedges, refprism_2fa_1fb_1ec_0v_splitfaces, 0, refprism_2fa_1fb_1ec_0v_newelstypes, refprism_2fa_1fb_1ec_0v_newels);

	//  HP_PRISM_2FA_1FB_2EB_0V  
	  public static int[][] refprism_2fa_1fb_2eb_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {0, 0, 0}
	  };

	 public static int[][] refprism_2fa_1fb_2eb_0v_splitfaces =
	 {
		 new int[] {2, 3, 5, 21},
		 new int[] {3, 2, 6, 22},
		 new int[] {3, 1, 6, 23},
		 new int[] {1, 3, 4, 24},
		 new int[] {5, 6, 2, 33},
		 new int[] {6, 5, 3, 34},
		 new int[] {6, 4, 3, 35},
		 new int[] {4, 1, 6, 36},
		 new int[] {4, 1, 5, 31},
		 new int[] {1, 2, 4, 19},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1fb_2eb_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1F_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_2fa_1fb_2eb_0v_newels =
	  {
		  new int[] {18, 23, 22, 30, 35, 34, 0, 0},
		  new int[] {24, 21, 22, 23, 36, 33, 34, 35},
		  new int[] {31, 29, 17, 19, 36, 33, 21, 24},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {7, 2, 9, 12, 19, 17, 21, 24},
		  new int[] {6, 43, 44, 30, 34, 35, 0, 0},
		  new int[] {44, 43, 42, 45, 35, 34, 33, 36},
		  new int[] {5, 40, 45, 42, 29, 31, 36, 33},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {40, 4, 45, 31, 28, 36, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_1fb_2eb_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1fb_2eb_0v_splitedges, refprism_2fa_1fb_2eb_0v_splitfaces, 0, refprism_2fa_1fb_2eb_0v_newelstypes, refprism_2fa_1fb_2eb_0v_newels);

	//  HP_PRISM_1FB_2EA_0V   ... quad face 1-2-4-5 with singular edges 1-4, 2-5
	  public static int[][] refprism_1fb_2ea_0v_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {1, 2, 9},
		  new int[] {2, 1, 10},
		  new int[] {4, 6, 11},
		  new int[] {5, 6, 12},
		  new int[] {4, 5, 13},
		  new int[] {5, 4, 14},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fb_2ea_0v_newelstypes = {HP_PRISM, HP_PRISM_1FB_1EA_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fb_2ea_0v_newels =
	  {
		  new int[] {7, 8, 3, 11, 12, 6, 0, 0},
		  new int[] {1, 9, 7, 4, 13, 11, 0, 0},
		  new int[] {13, 14, 10, 9, 11, 12, 8, 7},
		  new int[] {5, 14, 12, 2, 10, 8, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fb_2ea_0v = new HPRef_Struct(HP_PRISM, refprism_1fb_2ea_0v_splitedges, 0, 0, refprism_1fb_2ea_0v_newelstypes, refprism_1fb_2ea_0v_newels);

	//  HP_PRISM_1FB_2EB_0V   ... quad face 1-2-4-5 with singular edges 1-4, 3-6 
	  public static int[][] refprism_1fb_2eb_0v_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {2, 3, 9},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {1, 3, 12},
		  new int[] {4, 5, 40},
		  new int[] {5, 6, 42},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {4, 6, 45},
		  new int[] {0, 0, 0}
	  };
	public static HPREF_ELEMENT_TYPE[] refprism_1fb_2eb_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_PRISM_1FB_1EA_0V, HP_HEX_1F_0E_0V, HP_NONE};
	  public static int[][] refprism_1fb_2eb_0v_newels =
	  {
		  new int[] {3, 11, 10, 6, 44, 43, 0, 0},
		  new int[] {12, 9, 10, 11, 45, 42, 43, 44},
		  new int[] {1, 7, 12, 4, 40, 45, 0, 0},
		  new int[] {40, 5, 2, 7, 45, 42, 9, 12}
	  };
	  public static HPRef_Struct refprism_1fb_2eb_0v = new HPRef_Struct(HP_PRISM, refprism_1fb_2eb_0v_splitedges, 0, 0, refprism_1fb_2eb_0v_newelstypes, refprism_1fb_2eb_0v_newels);

	//  HP_PRISM_1FB_3E_0V   ... quad face 1-2-4-5 with singular edges 1-4, 3-6
	  public static int[][] refprism_1fb_3e_0v_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {2, 1, 8},
		  new int[] {2, 3, 9},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {1, 3, 12},
		  new int[] {4, 5, 40},
		  new int[] {5, 4, 41},
		  new int[] {5, 6, 42},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {4, 6, 45},
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refprism_1fb_3e_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_HEX, HP_PRISM_1FB_1EA_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fb_3e_0v_newels =
	  {
		  new int[] {3, 11, 10, 6, 44, 43, 0, 0},
		  new int[] {12, 9, 10, 11, 45, 42, 43, 44},
		  new int[] {1, 7, 12, 4, 40, 45, 0, 0},
		  new int[] {40, 41, 8, 7, 45, 42, 9, 12},
		  new int[] {5, 41, 42, 2, 8, 9, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fb_3e_0v = new HPRef_Struct(HP_PRISM, refprism_1fb_3e_0v_splitedges, 0, 0, refprism_1fb_3e_0v_newelstypes, refprism_1fb_3e_0v_newels);



	//  HP_PRISM_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3
	  public static int[][] refprism_2fb_0e_0v_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {1, 2, 9},
		  new int[] {3, 2, 10},
		  new int[] {4, 6, 11},
		  new int[] {5, 6, 12},
		  new int[] {4, 5, 13},
		  new int[] {6, 5, 14},
		  new int[] {0, 0, 0}
	  };
	 public static int[][] refprism_2fb_0e_0v_splitfaces =
	 {
		 new int[] {1, 2, 3, 15},
		 new int[] {4, 5, 6, 16},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_2fb_0e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_2fb_0e_0v_newels =
	  {
		  new int[] {15, 8, 10, 16, 12, 14, 0, 0},
		  new int[] {13, 5, 2, 9, 16, 12, 8, 15},
		  new int[] {11, 7, 3, 6, 16, 15, 10, 14},
		  new int[] {1, 9, 15, 4, 13, 16, 0, 0},
		  new int[] {4, 11, 16, 1, 7, 15, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fb_0e_0v = new HPRef_Struct(HP_PRISM, refprism_2fb_0e_0v_splitedges, refprism_2fb_0e_0v_splitfaces, 0, refprism_2fb_0e_0v_newelstypes, refprism_2fb_0e_0v_newels);

	//  HP_PRISM_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3 and sing edge 3-6
	  public static int[][] refprism_2fb_1ec_0v_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {1, 2, 9},
		  new int[] {3, 2, 10},
		  new int[] {4, 6, 11},
		  new int[] {5, 6, 12},
		  new int[] {4, 5, 13},
		  new int[] {6, 5, 14},
		  new int[] {3, 1, 17},
		  new int[] {6, 4, 18},
		  new int[] {0, 0, 0}
	  };
	 public static int[][] refprism_2fb_1ec_0v_splitfaces =
	 {
		 new int[] {1, 2, 3, 15},
		 new int[] {4, 5, 6, 16},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_2fb_1ec_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_2fb_1ec_0v_newels =
	  {
		  new int[] {15, 8, 10, 16, 12, 14, 0, 0},
		  new int[] {13, 5, 2, 9, 16, 12, 8, 15},
		  new int[] {11, 7, 17, 18, 16, 15, 10, 14},
		  new int[] {1, 9, 15, 4, 13, 16, 0, 0},
		  new int[] {4, 11, 16, 1, 7, 15, 0, 0},
		  new int[] {3, 17, 10, 6, 18, 14, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fb_1ec_0v = new HPRef_Struct(HP_PRISM, refprism_2fb_1ec_0v_splitedges, refprism_2fb_1ec_0v_splitfaces, 0, refprism_2fb_1ec_0v_newelstypes, refprism_2fb_1ec_0v_newels);



	//  HP_PRISM_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3 and 3 sing edges
	  public static int[][] refprism_2fb_3e_0v_splitedges =
	  {
		  new int[] {1, 3, 7},
		  new int[] {2, 3, 8},
		  new int[] {1, 2, 9},
		  new int[] {3, 2, 10},
		  new int[] {4, 6, 11},
		  new int[] {5, 6, 12},
		  new int[] {4, 5, 13},
		  new int[] {6, 5, 14},
		  new int[] {3, 1, 17},
		  new int[] {6, 4, 18},
		  new int[] {2, 1, 19},
		  new int[] {5, 4, 20},
		  new int[] {0, 0, 0}
	  };
	 public static int[][] refprism_2fb_3e_0v_splitfaces =
	 {
		 new int[] {1, 2, 3, 15},
		 new int[] {4, 5, 6, 16},
		 new int[] {0, 0, 0, 0}
	 };
	  public static HPREF_ELEMENT_TYPE[] refprism_2fb_3e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_2fb_3e_0v_newels =
	  {
		  new int[] {15, 8, 10, 16, 12, 14, 0, 0},
		  new int[] {13, 20, 19, 9, 16, 12, 8, 15},
		  new int[] {11, 7, 17, 18, 16, 15, 10, 14},
		  new int[] {1, 9, 15, 4, 13, 16, 0, 0},
		  new int[] {4, 11, 16, 1, 7, 15, 0, 0},
		  new int[] {3, 17, 10, 6, 18, 14, 0, 0},
		  new int[] {5, 20, 12, 2, 19, 8, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fb_3e_0v = new HPRef_Struct(HP_PRISM, refprism_2fb_3e_0v_splitedges, refprism_2fb_3e_0v_splitfaces, 0, refprism_2fb_3e_0v_newelstypes, refprism_2fb_3e_0v_newels);



	//  HP_PRISM_1FA_1FB_0E_0V   ... quad face 1-2-4-5 and trig face 1-2-3
	  public static int[][] refprism_1fa_1fb_0e_0v_splitedges =
	  {
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_1fb_0e_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {0, 0, 0, 0}
	  };

	public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_0e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_0e_0v_newels =
	  {
		  new int[] {24, 21, 18, 45, 42, 6, 0, 0},
		  new int[] {4, 5, 17, 16, 45, 42, 21, 24},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {1, 2, 9, 12, 16, 17, 21, 24}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_0e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_0e_0v_splitedges, refprism_1fa_1fb_0e_0v_splitfaces, 0, refprism_1fa_1fb_0e_0v_newelstypes, refprism_1fa_1fb_0e_0v_newels);

	/*
	//  HP_PRISM_1FA_1FB_1EC_0V   ... quad face 1-2-4-5 and trig face 1-2-3
	int refprism_1fa_1fb_1ec_0v_splitedges[][3] =
	    {
	      {1,4,16}, 
	      {2,5,17},
	      {3,6,18},
	      {2,3,9},
	      {1,3,12},
	      {5,6,42},
	      {4,6,45},
	      {6,5,43},
	      {6,4,44},
	      {3,2,10},
	      {3,1,11},
	      {0,0,0}
	    };
	  int refprism_1fa_1fb_1ec_0v_splitfaces[][4] = 
	    {
	      {2,3,5,21},
	      {1,3,4,24},
	      { 0, 0, 0, 0 }
	    };
	
	  HPREF_ELEMENT_TYPE refprism_1fa_1fb_1ec_0v_newelstypes[] =
	    {
	      HP_PRISM, 
	      HP_HEX_1F_0E_0V,
	      HP_PRISM_1FA_0E_0V, 
	      HP_HEX_1FA_1FB_0E_0V,
	      HP_PRISM_SINGEDGE,
	      HP_PRISM_1FA_1E_0V, 
	      HP_PRISM_
	      HP_NONE,
	    };
	  int refprism_1fa_1fb_0e_0v_newels[][8] =
	    {
	      { 24, 21, 18, 45, 42, 6 }, 
	      { 4, 5, 17, 16, 45, 42, 21, 24 },
	      { 12, 9, 3, 24, 21, 18 }, 
	      { 1, 2, 9, 12, 16, 17, 21, 24 } 
	    };
	  HPRef_Struct refprism_1fa_1fb_0e_0v =
	    {
	      HP_PRISM,
	      refprism_1fa_1fb_1ec_0v_splitedges, 
	
	      refprism_1fa_1fb_1ec_0v_splitfaces, 0,
	      refprism_1fa_1fb_1ec_0v_newelstypes, 
	      refprism_1fa_1fb_1ec_0v_newels
	    };
	
	
	*/




	//  HP_PRISM_2FA_1FB_0E_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
	  public static int[][] refprism_2fa_1fb_0e_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_2fa_1fb_0e_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {5, 6, 2, 33},
		  new int[] {4, 1, 6, 36},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1fb_0e_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PRISM, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_NONE};
	  public static int[][] refprism_2fa_1fb_0e_0v_newels =
	  {
		  new int[] {28, 29, 17, 16, 36, 33, 21, 24},
		  new int[] {24, 21, 18, 36, 33, 30, 0, 0},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {1, 2, 9, 12, 16, 17, 21, 24},
		  new int[] {6, 42, 45, 30, 33, 36, 0, 0},
		  new int[] {4, 5, 29, 28, 45, 42, 33, 36}
	  };
	  public static HPRef_Struct refprism_2fa_1fb_0e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1fb_0e_0v_splitedges, refprism_2fa_1fb_0e_0v_splitfaces, 0, refprism_2fa_1fb_0e_0v_newelstypes, refprism_2fa_1fb_0e_0v_newels);


	//  HP_PRISM_1FA_1FB_1EA_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
	  public static int[][] refprism_1fa_1fb_1ea_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_1fb_1ea_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {1, 2, 4, 19},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_1ea_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PRISM, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_1ea_0v_newels =
	  {
		  new int[] {40, 5, 17, 19, 45, 42, 21, 24},
		  new int[] {24, 21, 18, 45, 42, 6, 0, 0},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {7, 2, 9, 12, 19, 17, 21, 24},
		  new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_1ea_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_1ea_0v_splitedges, refprism_1fa_1fb_1ea_0v_splitfaces, 0, refprism_1fa_1fb_1ea_0v_newelstypes, refprism_1fa_1fb_1ea_0v_newels);

	//  HP_PRISM_2FA_1FB_1EA_0V   
	  public static int[][] refprism_2fa_1fb_1ea_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_2fa_1fb_1ea_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {1, 2, 4, 19},
		  new int[] {4, 1, 6, 36},
		  new int[] {4, 1, 5, 31},
		  new int[] {5, 6, 2, 33},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1fb_1ea_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_2fa_1fb_1ea_0v_newels =
	  {
		  new int[] {18, 24, 21, 30, 36, 33, 0, 0},
		  new int[] {31, 29, 17, 19, 36, 33, 21, 24},
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {3, 12, 9, 18, 24, 21, 0, 0},
		  new int[] {7, 2, 9, 12, 19, 17, 21, 24},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {6, 42, 45, 30, 33, 36, 0, 0},
		  new int[] {40, 5, 29, 31, 45, 42, 33, 36},
		  new int[] {40, 4, 45, 31, 28, 36, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_1fb_1ea_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1fb_1ea_0v_splitedges, refprism_2fa_1fb_1ea_0v_splitfaces, 0, refprism_2fa_1fb_1ea_0v_newelstypes, refprism_2fa_1fb_1ea_0v_newels);


	//  HP_PRISM_2FA_1FB_2EA_0V   
	  public static int[][] refprism_2fa_1fb_2ea_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_2fa_1fb_2ea_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {1, 2, 4, 19},
		  new int[] {4, 1, 6, 36},
		  new int[] {4, 1, 5, 31},
		  new int[] {5, 6, 2, 33},
		  new int[] {5, 4, 2, 32},
		  new int[] {2, 1, 5, 20},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1fb_2ea_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_2fa_1fb_2ea_0v_newels =
	  {
		  new int[] {18, 24, 21, 30, 36, 33, 0, 0},
		  new int[] {31, 32, 20, 19, 36, 33, 21, 24},
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {3, 12, 9, 18, 24, 21, 0, 0},
		  new int[] {7, 8, 9, 12, 19, 20, 21, 24},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {6, 42, 45, 30, 33, 36, 0, 0},
		  new int[] {40, 41, 32, 31, 45, 42, 33, 36},
		  new int[] {40, 4, 45, 31, 28, 36, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0},
		  new int[] {29, 32, 33, 17, 20, 21, 0, 0},
		  new int[] {5, 41, 42, 29, 32, 33, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_1fb_2ea_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1fb_2ea_0v_splitedges, refprism_2fa_1fb_2ea_0v_splitfaces, 0, refprism_2fa_1fb_2ea_0v_newelstypes, refprism_2fa_1fb_2ea_0v_newels);

	//  HP_PRISM_2FA_1FB_3E_0V   
	  public static int[][] refprism_2fa_1fb_3e_0v_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {2, 1, 8},
		  new int[] {2, 3, 9},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {4, 5, 40},
		  new int[] {5, 4, 41},
		  new int[] {5, 6, 42},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {4, 6, 45},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_2fa_1fb_3e_0v_splitfaces =
	  {
		  new int[] {1, 2, 4, 19},
		  new int[] {2, 1, 5, 20},
		  new int[] {2, 3, 5, 21},
		  new int[] {3, 2, 6, 22},
		  new int[] {3, 1, 6, 23},
		  new int[] {1, 3, 4, 24},
		  new int[] {4, 1, 5, 31},
		  new int[] {5, 4, 2, 32},
		  new int[] {5, 6, 2, 33},
		  new int[] {6, 5, 3, 34},
		  new int[] {6, 4, 3, 35},
		  new int[] {4, 1, 6, 36},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_1fb_3e_0v_newelstypes = {HP_HEX, HP_PRISM_SINGEDGE, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_HEX_1FA_1FB_0E_0V, HP_NONE};
	  public static int[][] refprism_2fa_1fb_3e_0v_newels =
	  {
		  new int[] {24, 21, 22, 23, 36, 33, 34, 35},
		  new int[] {18, 23, 22, 30, 35, 34, 0, 0},
		  new int[] {31, 32, 20, 19, 36, 33, 21, 24},
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {29, 32, 33, 17, 20, 21, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0},
		  new int[] {7, 8, 9, 12, 19, 20, 21, 24},
		  new int[] {44, 43, 42, 45, 35, 34, 33, 36},
		  new int[] {6, 43, 44, 30, 34, 35, 0, 0},
		  new int[] {40, 4, 45, 31, 28, 36, 0, 0},
		  new int[] {5, 41, 42, 29, 32, 33, 0, 0},
		  new int[] {40, 41, 32, 31, 45, 42, 33, 36}
	  };
	  public static HPRef_Struct refprism_2fa_1fb_3e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_1fb_3e_0v_splitedges, refprism_2fa_1fb_3e_0v_splitfaces, 0, refprism_2fa_1fb_3e_0v_newelstypes, refprism_2fa_1fb_3e_0v_newels);




	//  HP_PRISM_1FA_1FB_1EB_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
	  public static int[][] refprism_1fa_1fb_1eb_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_1fb_1eb_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {2, 1, 5, 20},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_1eb_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PRISM, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_1eb_0v_newels =
	  {
		  new int[] {4, 41, 20, 16, 45, 42, 21, 24},
		  new int[] {24, 21, 18, 45, 42, 6, 0, 0},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {1, 8, 9, 12, 16, 20, 21, 24},
		  new int[] {5, 41, 42, 17, 20, 21, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_1eb_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_1eb_0v_splitedges, refprism_1fa_1fb_1eb_0v_splitfaces, 0, refprism_1fa_1fb_1eb_0v_newelstypes, refprism_1fa_1fb_1eb_0v_newels);


	//  HP_PRISM_1FA_1FB_2EA_0V   ... quad face 1-2-4-5 and trig face 1-2-3 
	  public static int[][] refprism_1fa_1fb_2ea_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_1fb_2ea_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {2, 1, 5, 20},
		  new int[] {1, 2, 4, 19},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_2ea_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PRISM, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_2ea_0v_newels =
	  {
		  new int[] {40, 41, 20, 19, 45, 42, 21, 24},
		  new int[] {24, 21, 18, 45, 42, 6, 0, 0},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {7, 8, 9, 12, 19, 20, 21, 24},
		  new int[] {5, 41, 42, 17, 20, 21, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0},
		  new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_2ea_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_2ea_0v_splitedges, refprism_1fa_1fb_2ea_0v_splitfaces, 0, refprism_1fa_1fb_2ea_0v_newelstypes, refprism_1fa_1fb_2ea_0v_newels);


	//  HP_PRISM_1FA_1FB_3E_0V   
	  public static int[][] refprism_1fa_1fb_3e_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_1fb_3e_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {2, 1, 5, 20},
		  new int[] {1, 2, 4, 19},
		  new int[] {3, 2, 6, 22},
		  new int[] {3, 1, 6, 23},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_1fb_3e_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_HEX, HP_PRISM_SINGEDGE, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fa_1fb_3e_0v_newels =
	  {
		  new int[] {40, 41, 20, 19, 45, 42, 21, 24},
		  new int[] {24, 21, 22, 23, 45, 42, 43, 44},
		  new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		  new int[] {12, 9, 10, 11, 24, 21, 22, 23},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {7, 8, 9, 12, 19, 20, 21, 24},
		  new int[] {5, 41, 42, 17, 20, 21, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0},
		  new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_1fb_3e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_1fb_3e_0v_splitedges, refprism_1fa_1fb_3e_0v_splitfaces, 0, refprism_1fa_1fb_3e_0v_newelstypes, refprism_1fa_1fb_3e_0v_newels);








	//  HP_PRISM_2FA_0E_0V  singular trig faces
	  public static int[][] refprism_2fa_0e_0v_splitedges =
	  {
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {0, 0, 0}
	  };

	public static HPREF_ELEMENT_TYPE[] refprism_2fa_0e_0v_newelstypes = {HP_PRISM, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_NONE};
	  public static int[][] refprism_2fa_0e_0v_newels =
	  {
		  new int[] {16, 17, 18, 28, 29, 30, 0, 0},
		  new int[] {1, 2, 3, 16, 17, 18, 0, 0},
		  new int[] {4, 6, 5, 28, 30, 29, 0, 0}
	  };

	public static HPRef_Struct refprism_2fa_0e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_0e_0v_splitedges, 0, 0, refprism_2fa_0e_0v_newelstypes, refprism_2fa_0e_0v_newels);





	//  HP_PRISM_1FA_2FB    ... quad face 1-2-4-5 and quad face 1-4-6-3
	public static int[][] refprism_1fa_2fb_0e_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_0e_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_0e_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refprism_1fa_2fb_0e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_1fa_2fb_0e_0v_newels =
	  {
		  new int[] {25, 21, 22, 46, 42, 43, 0, 0},
		  new int[] {40, 5, 17, 19, 46, 42, 21, 25},
		  new int[] {24, 18, 6, 45, 25, 22, 43, 46},
		  new int[] {16, 19, 25, 4, 40, 46, 0, 0},
		  new int[] {4, 45, 46, 16, 24, 25, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 2, 9, 13, 19, 17, 21, 25},
		  new int[] {3, 12, 13, 10, 18, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_2fb_0e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_2fb_0e_0v_splitedges, refprism_1fa_2fb_0e_0v_splitfaces, refprism_1fa_2fb_0e_0v_splitelement, refprism_1fa_2fb_0e_0v_newelstypes, refprism_1fa_2fb_0e_0v_newels);

	//  HP_PRISM_1FA_2FB_1EC    ... quad face 1-2-4-5 and quad face 1-4-6-3
	public static int[][] refprism_1fa_2fb_1ec_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_1ec_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_1ec_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refprism_1fa_2fb_1ec_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_1fa_2fb_1ec_0v_newels =
	  {
		  new int[] {25, 21, 22, 46, 42, 43, 0, 0},
		  new int[] {40, 5, 17, 19, 46, 42, 21, 25},
		  new int[] {24, 23, 44, 45, 25, 22, 43, 46},
		  new int[] {16, 19, 25, 4, 40, 46, 0, 0},
		  new int[] {4, 45, 46, 16, 24, 25, 0, 0},
		  new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 2, 9, 13, 19, 17, 21, 25},
		  new int[] {11, 12, 13, 10, 23, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_2fb_1ec_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_2fb_1ec_0v_splitedges, refprism_1fa_2fb_1ec_0v_splitfaces, refprism_1fa_2fb_1ec_0v_splitelement, refprism_1fa_2fb_1ec_0v_newelstypes, refprism_1fa_2fb_1ec_0v_newels);


	//  HP_PRISM_1FA_2FB_3E    ... quad face 1-2-4-5 and quad face 1-4-6-3
	public static int[][] refprism_1fa_2fb_3e_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_3e_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_3e_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refprism_1fa_2fb_3e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_1fa_2fb_3e_0v_newels =
	  {
		  new int[] {25, 21, 22, 46, 42, 43, 0, 0},
		  new int[] {40, 41, 20, 19, 46, 42, 21, 25},
		  new int[] {24, 23, 44, 45, 25, 22, 43, 46},
		  new int[] {16, 19, 25, 4, 40, 46, 0, 0},
		  new int[] {4, 45, 46, 16, 24, 25, 0, 0},
		  new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		  new int[] {5, 41, 42, 17, 20, 21, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 8, 9, 13, 19, 20, 21, 25},
		  new int[] {11, 12, 13, 10, 23, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_2fb_3e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_2fb_3e_0v_splitedges, refprism_1fa_2fb_3e_0v_splitfaces, refprism_1fa_2fb_3e_0v_splitelement, refprism_1fa_2fb_3e_0v_newelstypes, refprism_1fa_2fb_3e_0v_newels);









	//  HP_PRISM_1FA_2FB_1eb    ... quad face 1-2-4-5 and quad face 1-4-6-3
	public static int[][] refprism_1fa_2fb_1eb_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_1eb_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_1fa_2fb_1eb_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {0, 0, 0, 0, 0}
	};


	public static HPREF_ELEMENT_TYPE[] refprism_1fa_2fb_1eb_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};

	  public static int[][] refprism_1fa_2fb_1eb_0v_newels =
	  {
		  new int[] {25, 21, 22, 46, 42, 43, 0, 0},
		  new int[] {40, 41, 20, 19, 46, 42, 21, 25},
		  new int[] {24, 18, 6, 45, 25, 22, 43, 46},
		  new int[] {16, 19, 25, 4, 40, 46, 0, 0},
		  new int[] {4, 45, 46, 16, 24, 25, 0, 0},
		  new int[] {5, 41, 42, 17, 20, 21, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 8, 9, 13, 19, 20, 21, 25},
		  new int[] {3, 12, 13, 10, 18, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_2fb_1eb_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_2fb_1eb_0v_splitedges, refprism_1fa_2fb_1eb_0v_splitfaces, refprism_1fa_2fb_1eb_0v_splitelement, refprism_1fa_2fb_1eb_0v_newelstypes, refprism_1fa_2fb_1eb_0v_newels);






	//  HP_PRISM_2FA_2FB 
	public static int[][] refprism_2fa_2fb_0e_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {4, 6, 45},
		new int[] {4, 1, 28},
		new int[] {5, 2, 29},
		new int[] {6, 3, 30},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_2fa_2fb_0e_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {4, 1, 5, 31},
		new int[] {5, 6, 2, 33},
		new int[] {6, 5, 3, 34},
		new int[] {4, 1, 6, 36},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_2fa_2fb_0e_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {4, 1, 6, 5, 37},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refprism_2fa_2fb_0e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_2fa_2fb_0e_0v_newels =
	  {
		  new int[] {25, 21, 22, 37, 33, 34, 0, 0},
		  new int[] {31, 29, 17, 19, 37, 33, 21, 25},
		  new int[] {36, 24, 18, 30, 37, 25, 22, 34},
		  new int[] {16, 19, 25, 28, 31, 37, 0, 0},
		  new int[] {28, 36, 37, 16, 24, 25, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 2, 9, 13, 19, 17, 21, 25},
		  new int[] {3, 12, 13, 10, 18, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {46, 43, 42, 37, 34, 33, 0, 0},
		  new int[] {40, 5, 29, 31, 46, 42, 33, 37},
		  new int[] {6, 45, 36, 30, 43, 46, 37, 34},
		  new int[] {40, 4, 46, 31, 28, 37, 0, 0},
		  new int[] {4, 45, 46, 28, 36, 37, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_2fb_0e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_2fb_0e_0v_splitedges, refprism_2fa_2fb_0e_0v_splitfaces, refprism_2fa_2fb_0e_0v_splitelement, refprism_2fa_2fb_0e_0v_newelstypes, refprism_2fa_2fb_0e_0v_newels);


	//  HP_PRISM_2FA_2FB_1EC 
	public static int[][] refprism_2fa_2fb_1ec_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 1, 28},
		new int[] {5, 2, 29},
		new int[] {6, 3, 30},
		new int[] {4, 5, 40},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_2fa_2fb_1ec_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {4, 1, 5, 31},
		new int[] {5, 6, 2, 33},
		new int[] {6, 5, 3, 34},
		new int[] {6, 4, 3, 35},
		new int[] {4, 1, 6, 36},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_2fa_2fb_1ec_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {4, 1, 6, 5, 37},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refprism_2fa_2fb_1ec_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_2fa_2fb_1ec_0v_newels =
	  {
		  new int[] {25, 21, 22, 37, 33, 34, 0, 0},
		  new int[] {31, 29, 17, 19, 37, 33, 21, 25},
		  new int[] {36, 24, 23, 35, 37, 25, 22, 34},
		  new int[] {16, 19, 25, 28, 31, 37, 0, 0},
		  new int[] {28, 36, 37, 16, 24, 25, 0, 0},
		  new int[] {18, 23, 22, 30, 35, 34, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 2, 9, 13, 19, 17, 21, 25},
		  new int[] {11, 12, 13, 10, 23, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {46, 43, 42, 37, 34, 33, 0, 0},
		  new int[] {40, 5, 29, 31, 46, 42, 33, 37},
		  new int[] {44, 45, 36, 35, 43, 46, 37, 34},
		  new int[] {40, 4, 46, 31, 28, 37, 0, 0},
		  new int[] {4, 45, 46, 28, 36, 37, 0, 0},
		  new int[] {44, 6, 43, 35, 30, 34, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_2fb_1ec_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_2fb_1ec_0v_splitedges, refprism_2fa_2fb_1ec_0v_splitfaces, refprism_2fa_2fb_1ec_0v_splitelement, refprism_2fa_2fb_1ec_0v_newelstypes, refprism_2fa_2fb_1ec_0v_newels);



	//  HP_PRISM_2FA_2FB_3E 
	public static int[][] refprism_2fa_2fb_3e_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 1, 28},
		new int[] {5, 2, 29},
		new int[] {6, 3, 30},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_2fa_2fb_3e_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {4, 1, 5, 31},
		new int[] {5, 4, 2, 32},
		new int[] {5, 6, 2, 33},
		new int[] {6, 5, 3, 34},
		new int[] {6, 4, 3, 35},
		new int[] {4, 1, 6, 36},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] refprism_2fa_2fb_3e_0v_splitelement =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {4, 1, 6, 5, 37},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refprism_2fa_2fb_3e_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_2fa_2fb_3e_0v_newels =
	  {
		  new int[] {25, 21, 22, 37, 33, 34, 0, 0},
		  new int[] {31, 32, 20, 19, 37, 33, 21, 25},
		  new int[] {36, 24, 23, 35, 37, 25, 22, 34},
		  new int[] {16, 19, 25, 28, 31, 37, 0, 0},
		  new int[] {28, 36, 37, 16, 24, 25, 0, 0},
		  new int[] {18, 23, 22, 30, 35, 34, 0, 0},
		  new int[] {29, 32, 33, 17, 20, 21, 0, 0},
		  new int[] {13, 9, 10, 25, 21, 22, 0, 0},
		  new int[] {7, 8, 9, 13, 19, 20, 21, 25},
		  new int[] {11, 12, 13, 10, 23, 24, 25, 22},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {8, 2, 9, 20, 17, 21, 0, 0},
		  new int[] {46, 43, 42, 37, 34, 33, 0, 0},
		  new int[] {40, 41, 32, 31, 46, 42, 33, 37},
		  new int[] {44, 45, 36, 35, 43, 46, 37, 34},
		  new int[] {40, 4, 46, 31, 28, 37, 0, 0},
		  new int[] {4, 45, 46, 28, 36, 37, 0, 0},
		  new int[] {44, 6, 43, 35, 30, 34, 0, 0},
		  new int[] {5, 41, 42, 29, 32, 33, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_2fb_3e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_2fb_3e_0v_splitedges, refprism_2fa_2fb_3e_0v_splitfaces, refprism_2fa_2fb_3e_0v_splitelement, refprism_2fa_2fb_3e_0v_newelstypes, refprism_2fa_2fb_3e_0v_newels);




	//  HP_PRISM_1FA_2E_0V  
	  public static int[][] refprism_1fa_2e_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_1fa_2e_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {2, 1, 5, 20},
		  new int[] {1, 2, 4, 19},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_2e_0v_newelstypes = {HP_HEX, HP_PRISM, HP_PRISM_1FA_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_SINGEDGE, HP_PRISM_1FA_1E_0V, HP_PRISM_SINGEDGE, HP_PRISM_1FA_1E_0V, HP_NONE};
	  public static int[][] refprism_1fa_2e_0v_newels =
	  {
		  new int[] {40, 41, 20, 19, 45, 42, 21, 24},
		  new int[] {24, 21, 18, 45, 42, 6, 0, 0},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {9, 12, 7, 8, 21, 24, 19, 20},
		  new int[] {17, 21, 20, 5, 42, 41, 0, 0},
		  new int[] {2, 9, 8, 17, 21, 20, 0, 0},
		  new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_2e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_2e_0v_splitedges, refprism_1fa_2e_0v_splitfaces, 0, refprism_1fa_2e_0v_newelstypes, refprism_1fa_2e_0v_newels);

	//  HP_PRISM_2FA_2E_0V   
	  public static int[][] refprism_2fa_2e_0v_splitedges =
	  {
		  new int[] {2, 3, 9},
		  new int[] {1, 3, 12},
		  new int[] {1, 4, 16},
		  new int[] {2, 5, 17},
		  new int[] {3, 6, 18},
		  new int[] {5, 6, 42},
		  new int[] {4, 6, 45},
		  new int[] {4, 1, 28},
		  new int[] {5, 2, 29},
		  new int[] {6, 3, 30},
		  new int[] {4, 5, 40},
		  new int[] {1, 2, 7},
		  new int[] {5, 4, 41},
		  new int[] {2, 1, 8},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_2fa_2e_0v_splitfaces =
	  {
		  new int[] {2, 3, 5, 21},
		  new int[] {1, 3, 4, 24},
		  new int[] {1, 2, 4, 19},
		  new int[] {4, 1, 6, 36},
		  new int[] {4, 1, 5, 31},
		  new int[] {5, 6, 2, 33},
		  new int[] {5, 4, 2, 32},
		  new int[] {2, 1, 5, 20},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_2e_0v_newelstypes = {HP_PRISM, HP_HEX, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_1FA_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_NONE};
	  public static int[][] refprism_2fa_2e_0v_newels =
	  {
		  new int[] {24, 21, 18, 36, 33, 30, 0, 0},
		  new int[] {19, 20, 21, 24, 31, 32, 33, 36},
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {17, 21, 20, 29, 33, 32, 0, 0},
		  new int[] {12, 9, 3, 24, 21, 18, 0, 0},
		  new int[] {7, 8, 9, 12, 19, 20, 21, 24},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {2, 9, 8, 17, 21, 20, 0, 0},
		  new int[] {45, 6, 42, 36, 30, 33, 0, 0},
		  new int[] {40, 45, 42, 41, 31, 36, 33, 32},
		  new int[] {4, 45, 40, 28, 36, 31, 0, 0},
		  new int[] {5, 41, 42, 29, 32, 33, 0, 0}
	  };
	  public static HPRef_Struct refprism_2fa_2e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_2e_0v_splitedges, refprism_2fa_2e_0v_splitfaces, 0, refprism_2fa_2e_0v_newelstypes, refprism_2fa_2e_0v_newels);



	//  HP_PRISM_3E_0V   
	  public static int[][] refprism_3e_0v_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {2, 1, 8},
		  new int[] {2, 3, 9},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {1, 3, 12},
		  new int[] {4, 5, 40},
		  new int[] {5, 4, 41},
		  new int[] {5, 6, 42},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {4, 6, 45},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_3e_0v_splitfaces =
	  {
		  new int[] {1, 2, 3, 13},
		  new int[] {2, 3, 1, 14},
		  new int[] {3, 1, 2, 15},
		  new int[] {4, 5, 6, 46},
		  new int[] {5, 4, 6, 47},
		  new int[] {6, 4, 5, 48},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_3e_0v_newelstypes = {HP_PRISM, HP_HEX, HP_HEX, HP_HEX, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_NONE};
	  public static int[][] refprism_3e_0v_newels =
	  {
		  new int[] {13, 14, 15, 46, 47, 48, 0, 0},
		  new int[] {7, 8, 14, 13, 40, 41, 47, 46},
		  new int[] {15, 14, 9, 10, 48, 47, 42, 43},
		  new int[] {12, 13, 15, 11, 45, 46, 48, 44},
		  new int[] {14, 8, 9, 47, 41, 42, 0, 0},
		  new int[] {11, 15, 10, 44, 48, 43, 0, 0},
		  new int[] {7, 13, 12, 40, 46, 45, 0, 0},
		  new int[] {1, 7, 12, 4, 40, 45, 0, 0},
		  new int[] {2, 9, 8, 5, 42, 41, 0, 0},
		  new int[] {3, 11, 10, 6, 44, 43, 0, 0}
	  };
	  public static HPRef_Struct refprism_3e_0v = new HPRef_Struct(HP_PRISM, refprism_3e_0v_splitedges, refprism_3e_0v_splitfaces, 0, refprism_3e_0v_newelstypes, refprism_3e_0v_newels);


	//  HP_PRISM_3E_0V   
	public static int[][] refprism_1fa_3e_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_1fa_3e_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {2, 3, 1, 14},
		new int[] {3, 1, 2, 15},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {5, 4, 6, 47},
		new int[] {6, 4, 5, 48},
		new int[] {0, 0, 0, 0}
	};

	public static int[][] refprism_1fa_3e_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {2, 1, 3, 5, 26},
		new int[] {3, 1, 2, 6, 27},
		new int[] {0, 0, 0, 0, 0}
	};

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_3e_0v_newelstypes = {HP_PRISM, HP_HEX, HP_HEX, HP_HEX, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_1FA_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_NONE};
	public static int[][] refprism_1fa_3e_0v_newels =
	{
		new int[] {25, 26, 27, 46, 47, 48, 0, 0},
		new int[] {19, 20, 26, 25, 40, 41, 47, 46},
		new int[] {27, 26, 21, 22, 48, 47, 42, 43},
		new int[] {23, 24, 25, 27, 44, 45, 46, 48},
		new int[] {19, 25, 24, 40, 46, 45, 0, 0},
		new int[] {26, 20, 21, 47, 41, 42, 0, 0},
		new int[] {23, 27, 22, 44, 48, 43, 0, 0},
		new int[] {16, 19, 24, 4, 40, 45, 0, 0},
		new int[] {17, 21, 20, 5, 42, 41, 0, 0},
		new int[] {18, 23, 22, 6, 44, 43, 0, 0},
		new int[] {13, 14, 15, 25, 26, 27, 0, 0},
		new int[] {7, 8, 14, 13, 19, 20, 26, 25},
		new int[] {15, 14, 9, 10, 27, 26, 21, 22},
		new int[] {12, 13, 15, 11, 24, 25, 27, 23},
		new int[] {14, 8, 9, 26, 20, 21, 0, 0},
		new int[] {11, 15, 10, 23, 27, 22, 0, 0},
		new int[] {7, 13, 12, 19, 25, 24, 0, 0},
		new int[] {2, 9, 8, 17, 21, 20, 0, 0},
		new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		new int[] {1, 7, 12, 16, 19, 24, 0, 0}
	};
	  public static HPRef_Struct refprism_1fa_3e_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_3e_0v_splitedges, refprism_1fa_3e_0v_splitfaces, refprism_1fa_3e_0v_splitelements, refprism_1fa_3e_0v_newelstypes, refprism_1fa_3e_0v_newels);



	//  HP_PRISM_2FA_3E_0V   
	public static int[][] refprism_2fa_3e_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 1, 28},
		new int[] {5, 2, 29},
		new int[] {6, 3, 30},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_2fa_3e_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {2, 3, 1, 14},
		new int[] {3, 1, 2, 15},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 1, 5, 31},
		new int[] {5, 4, 2, 32},
		new int[] {5, 6, 2, 33},
		new int[] {6, 5, 3, 34},
		new int[] {6, 4, 3, 35},
		new int[] {4, 1, 6, 36},
		new int[] {4, 5, 6, 46},
		new int[] {5, 4, 6, 47},
		new int[] {6, 4, 5, 48},
		new int[] {0, 0, 0, 0}
	};

	public static int[][] refprism_2fa_3e_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {2, 1, 3, 5, 26},
		new int[] {3, 1, 2, 6, 27},
		new int[] {4, 1, 6, 5, 37},
		new int[] {5, 2, 4, 6, 38},
		new int[] {6, 4, 5, 3, 39},
		new int[] {0, 0, 0, 0, 0}
	};

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_3e_0v_newelstypes = {HP_PRISM, HP_HEX, HP_HEX, HP_HEX, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_1FA_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_PRISM_1FA_1E_0V, HP_NONE};

	  public static int[][] refprism_2fa_3e_0v_newels =
	  {
		  new int[] {25, 26, 27, 37, 38, 39, 0, 0},
		  new int[] {19, 20, 26, 25, 31, 32, 38, 37},
		  new int[] {27, 26, 21, 22, 39, 38, 33, 34},
		  new int[] {23, 24, 25, 27, 35, 36, 37, 39},
		  new int[] {19, 25, 24, 31, 37, 36, 0, 0},
		  new int[] {26, 20, 21, 38, 32, 33, 0, 0},
		  new int[] {23, 27, 22, 35, 39, 34, 0, 0},
		  new int[] {16, 19, 24, 28, 31, 36, 0, 0},
		  new int[] {17, 21, 20, 29, 33, 32, 0, 0},
		  new int[] {18, 23, 22, 30, 35, 34, 0, 0},
		  new int[] {13, 14, 15, 25, 26, 27, 0, 0},
		  new int[] {7, 8, 14, 13, 19, 20, 26, 25},
		  new int[] {15, 14, 9, 10, 27, 26, 21, 22},
		  new int[] {12, 13, 15, 11, 24, 25, 27, 23},
		  new int[] {14, 8, 9, 26, 20, 21, 0, 0},
		  new int[] {11, 15, 10, 23, 27, 22, 0, 0},
		  new int[] {7, 13, 12, 19, 25, 24, 0, 0},
		  new int[] {2, 9, 8, 17, 21, 20, 0, 0},
		  new int[] {3, 11, 10, 18, 23, 22, 0, 0},
		  new int[] {1, 7, 12, 16, 19, 24, 0, 0},
		  new int[] {48, 47, 46, 39, 38, 37, 0, 0},
		  new int[] {48, 43, 42, 47, 39, 34, 33, 38},
		  new int[] {45, 44, 48, 46, 36, 35, 39, 37},
		  new int[] {46, 47, 41, 40, 37, 38, 32, 31},
		  new int[] {47, 42, 41, 38, 33, 32, 0, 0},
		  new int[] {45, 46, 40, 36, 37, 31, 0, 0},
		  new int[] {44, 43, 48, 35, 34, 39, 0, 0},
		  new int[] {6, 43, 44, 30, 34, 35, 0, 0},
		  new int[] {5, 41, 42, 29, 32, 33, 0, 0},
		  new int[] {4, 45, 40, 28, 36, 31, 0, 0}
	  };

	public static HPRef_Struct refprism_2fa_3e_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_3e_0v_splitedges, refprism_2fa_3e_0v_splitfaces, refprism_2fa_3e_0v_splitelements, refprism_2fa_3e_0v_newelstypes, refprism_2fa_3e_0v_newels);



	//  HP_PRISM_3FB_0V   
	  public static int[][] refprism_3fb_0v_splitedges =
	  {
		  new int[] {1, 2, 7},
		  new int[] {2, 1, 8},
		  new int[] {2, 3, 9},
		  new int[] {3, 2, 10},
		  new int[] {3, 1, 11},
		  new int[] {1, 3, 12},
		  new int[] {4, 5, 40},
		  new int[] {5, 4, 41},
		  new int[] {5, 6, 42},
		  new int[] {6, 5, 43},
		  new int[] {6, 4, 44},
		  new int[] {4, 6, 45},
		  new int[] {0, 0, 0}
	  };
	  public static int[][] refprism_3fb_0v_splitfaces =
	  {
		  new int[] {1, 2, 3, 13},
		  new int[] {2, 3, 1, 14},
		  new int[] {3, 1, 2, 15},
		  new int[] {4, 5, 6, 46},
		  new int[] {5, 4, 6, 47},
		  new int[] {6, 4, 5, 48},
		  new int[] {0, 0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refprism_3fb_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_NONE};
	  public static int[][] refprism_3fb_0v_newels =
	  {
		  new int[] {13, 14, 15, 46, 47, 48, 0, 0},
		  new int[] {8, 7, 40, 41, 14, 13, 46, 47},
		  new int[] {10, 9, 42, 43, 15, 14, 47, 48},
		  new int[] {44, 45, 12, 11, 48, 46, 13, 15},
		  new int[] {1, 7, 13, 4, 40, 46, 0, 0},
		  new int[] {4, 45, 46, 1, 12, 13, 0, 0},
		  new int[] {2, 9, 14, 5, 42, 47, 0, 0},
		  new int[] {5, 41, 47, 2, 8, 14, 0, 0},
		  new int[] {3, 11, 15, 6, 44, 48, 0, 0},
		  new int[] {6, 43, 48, 3, 10, 15, 0, 0}
	  };
	  public static HPRef_Struct refprism_3fb_0v = new HPRef_Struct(HP_PRISM, refprism_3fb_0v_splitedges, refprism_3fb_0v_splitfaces, 0, refprism_3fb_0v_newelstypes, refprism_3fb_0v_newels);


	//  HP_PRISM_3FB_0V   
	public static int[][] refprism_1fa_3fb_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_1fa_3fb_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {2, 3, 1, 14},
		new int[] {3, 1, 2, 15},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 5, 6, 46},
		new int[] {5, 4, 6, 47},
		new int[] {6, 4, 5, 48},
		new int[] {0, 0, 0, 0}
	};

	public static int[][] refprism_1fa_3fb_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {2, 1, 3, 5, 26},
		new int[] {3, 1, 2, 6, 27},
		new int[] {0, 0, 0, 0, 0}
	};

	  public static HPREF_ELEMENT_TYPE[] refprism_1fa_3fb_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_1fa_3fb_0v_newels =
	  {
		  new int[] {25, 26, 27, 46, 47, 48, 0, 0},
		  new int[] {19, 40, 41, 20, 25, 46, 47, 26},
		  new int[] {22, 21, 42, 43, 27, 26, 47, 48},
		  new int[] {24, 23, 44, 45, 25, 27, 48, 46},
		  new int[] {16, 19, 25, 4, 40, 46, 0, 0},
		  new int[] {4, 45, 46, 16, 24, 25, 0, 0},
		  new int[] {17, 21, 26, 5, 42, 47, 0, 0},
		  new int[] {5, 41, 47, 17, 20, 26, 0, 0},
		  new int[] {18, 23, 27, 6, 44, 48, 0, 0},
		  new int[] {6, 43, 48, 18, 22, 27, 0, 0},
		  new int[] {13, 14, 15, 25, 26, 27, 0, 0},
		  new int[] {7, 8, 14, 13, 19, 20, 26, 25},
		  new int[] {9, 10, 15, 14, 21, 22, 27, 26},
		  new int[] {11, 12, 13, 15, 23, 24, 25, 27},
		  new int[] {2, 9, 14, 17, 21, 26, 0, 0},
		  new int[] {8, 2, 14, 20, 17, 26, 0, 0},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {3, 11, 15, 18, 23, 27, 0, 0},
		  new int[] {10, 3, 15, 22, 18, 27, 0, 0}
	  };
	  public static HPRef_Struct refprism_1fa_3fb_0v = new HPRef_Struct(HP_PRISM, refprism_1fa_3fb_0v_splitedges, refprism_1fa_3fb_0v_splitfaces, refprism_1fa_3fb_0v_splitelements, refprism_1fa_3fb_0v_newelstypes, refprism_1fa_3fb_0v_newels);



	//  HP_PRISM_2FA_3E_0V   
	public static int[][] refprism_2fa_3fb_0v_splitedges =
	{
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {3, 2, 10},
		new int[] {3, 1, 11},
		new int[] {1, 3, 12},
		new int[] {1, 4, 16},
		new int[] {2, 5, 17},
		new int[] {3, 6, 18},
		new int[] {4, 1, 28},
		new int[] {5, 2, 29},
		new int[] {6, 3, 30},
		new int[] {4, 5, 40},
		new int[] {5, 4, 41},
		new int[] {5, 6, 42},
		new int[] {6, 5, 43},
		new int[] {6, 4, 44},
		new int[] {4, 6, 45},
		new int[] {0, 0, 0}
	};
	public static int[][] refprism_2fa_3fb_0v_splitfaces =
	{
		new int[] {1, 2, 3, 13},
		new int[] {2, 3, 1, 14},
		new int[] {3, 1, 2, 15},
		new int[] {1, 2, 4, 19},
		new int[] {2, 1, 5, 20},
		new int[] {2, 3, 5, 21},
		new int[] {3, 2, 6, 22},
		new int[] {3, 1, 6, 23},
		new int[] {1, 3, 4, 24},
		new int[] {4, 1, 5, 31},
		new int[] {5, 4, 2, 32},
		new int[] {5, 6, 2, 33},
		new int[] {6, 5, 3, 34},
		new int[] {6, 4, 3, 35},
		new int[] {4, 1, 6, 36},
		new int[] {4, 5, 6, 46},
		new int[] {5, 4, 6, 47},
		new int[] {6, 4, 5, 48},
		new int[] {0, 0, 0, 0}
	};

	public static int[][] refprism_2fa_3fb_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 25},
		new int[] {2, 1, 3, 5, 26},
		new int[] {3, 1, 2, 6, 27},
		new int[] {4, 1, 6, 5, 37},
		new int[] {5, 2, 4, 6, 38},
		new int[] {6, 4, 5, 3, 39},
		new int[] {0, 0, 0, 0, 0}
	};

	  public static HPREF_ELEMENT_TYPE[] refprism_2fa_3fb_0v_newelstypes = {HP_PRISM, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_HEX_1F_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_HEX_1FA_1FB_0E_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_PRISM_1FA_1FB_1EA_0V, HP_PRISM_1FA_1FB_1EB_0V, HP_NONE};
	  public static int[][] refprism_2fa_3fb_0v_newels =
	  {
		  new int[] {25, 26, 27, 37, 38, 39, 0, 0},
		  new int[] {19, 31, 32, 20, 25, 37, 38, 26},
		  new int[] {33, 34, 22, 21, 38, 39, 27, 26},
		  new int[] {35, 36, 24, 23, 39, 37, 25, 27},
		  new int[] {16, 19, 25, 28, 31, 37, 0, 0},
		  new int[] {28, 36, 37, 16, 24, 25, 0, 0},
		  new int[] {17, 21, 26, 29, 33, 38, 0, 0},
		  new int[] {29, 32, 38, 17, 20, 26, 0, 0},
		  new int[] {18, 23, 27, 30, 35, 39, 0, 0},
		  new int[] {30, 34, 39, 18, 22, 27, 0, 0},
		  new int[] {13, 14, 15, 25, 26, 27, 0, 0},
		  new int[] {7, 8, 14, 13, 19, 20, 26, 25},
		  new int[] {9, 10, 15, 14, 21, 22, 27, 26},
		  new int[] {11, 12, 13, 15, 23, 24, 25, 27},
		  new int[] {2, 9, 14, 17, 21, 26, 0, 0},
		  new int[] {8, 2, 14, 20, 17, 26, 0, 0},
		  new int[] {1, 7, 13, 16, 19, 25, 0, 0},
		  new int[] {12, 1, 13, 24, 16, 25, 0, 0},
		  new int[] {3, 11, 15, 18, 23, 27, 0, 0},
		  new int[] {10, 3, 15, 22, 18, 27, 0, 0},
		  new int[] {48, 47, 46, 39, 38, 37, 0, 0},
		  new int[] {44, 45, 36, 35, 48, 46, 37, 39},
		  new int[] {40, 41, 32, 31, 46, 47, 38, 37},
		  new int[] {42, 43, 34, 33, 47, 48, 39, 38},
		  new int[] {6, 43, 48, 30, 34, 39, 0, 0},
		  new int[] {44, 6, 48, 35, 30, 39, 0, 0},
		  new int[] {4, 45, 46, 28, 36, 37, 0, 0},
		  new int[] {40, 4, 46, 31, 28, 37, 0, 0},
		  new int[] {5, 41, 47, 29, 32, 38, 0, 0},
		  new int[] {42, 5, 47, 33, 29, 38, 0, 0}
	  };

	public static HPRef_Struct refprism_2fa_3fb_0v = new HPRef_Struct(HP_PRISM, refprism_2fa_3fb_0v_splitedges, refprism_2fa_3fb_0v_splitfaces, refprism_2fa_3fb_0v_splitelements, refprism_2fa_3fb_0v_newelstypes, refprism_2fa_3fb_0v_newels);


	/* 
	
	
	//  HP_PRISM_3E_4EH
	int refprism_3e_4eh_splitedges[][3] =
	    {
	      { 1, 2, 7},
	      { 2, 1, 8},
	      { 2, 3, 9},
	      { 3, 2, 10},
	      { 3, 1, 11},
	      { 1, 3, 12},
	      { 4, 5, 40},
	      { 5, 4, 41},
	      { 5, 6, 42},
	      { 6, 5, 43},
	      { 6, 4, 44},
	      { 4, 6, 45},
	      { 0, 0, 0},
	
	    };
	int refprism_3e_4eh_splitfaces[][4] = 
	    {
	      {3,1,2,15},
	      {6,4,5,48}, 
	      {0,0,0,0}, 
	    };
	
	HPREF_ELEMENT_TYPE refprism_2fa_3fb_0v_newelstypes[] =
	  {
	    HP_PRISM, 
	    HP_HEX_2EH_0V,
	    HP_HEX_2EH_0V,
	    HP_TET_2E,
	    HP_TET_2E,
	    HP_PRISM_1E_2EH_0V, 
	    HP_PRISM_1E_2EH_0V, 
	    HP_NONE
	    };
	  int refprism_2fa_3fb_0v_newels[][8] =
	    {
	      {15, 7, 8, 48, 40, 41 }, 
	      
	    };
	
	HPRef_Struct refprism_2fa_3fb_0v =
	  {
	    HP_PRISM,
	    refprism_2fa_3fb_0v_splitedges, 
	    refprism_2fa_3fb_0v_splitfaces, 
	    refprism_2fa_3fb_0v_splitelements, 
	    refprism_2fa_3fb_0v_newelstypes, 
	    refprism_2fa_3fb_0v_newels
	  };
	*/ 

	/*
	//  HP_PRISM_2FA_3E_0V   
	int refprism_3e_4_0v_splitedges[][3] =
	    {
	      { 1, 2, 7},
	      { 2, 1, 8},
	      { 2, 3, 9},
	      { 3, 2, 10},
	      { 3, 1, 11},
	      { 1, 3, 12},
	      { 1, 4, 16}, 
	      { 2, 5, 17},
	      { 3, 6, 18},
	      { 4, 1, 28},
	      { 5, 2, 29},
	      { 6, 3, 30},
	      { 4, 5, 40},
	      { 5, 4, 41},
	      { 5, 6, 42},
	      { 6, 5, 43},
	      { 6, 4, 44},
	      { 4, 6, 45},
	      { 0, 0, 0}, 
	    };
	int refprism_2fa_3e_0v_splitfaces[][4] = 
	    {
	      {1,2,3,13},
	      {2,3,1,14},
	      {3,1,2,15},
	      {1,2,4,19},
	      {2,1,5,20},
	      {2,3,5,21},
	      {3,2,6,22},
	      {3,1,6,23},
	      {1,3,4,24},
	      {4,1,5,31},
	      {5,4,2,32},
	      {5,6,2,33},
	      {6,5,3,34},
	      {6,4,3,35},
	      {4,1,6,36},
	      {4,5,6,46},
	      {5,4,6,47},
	      {6,4,5,48}, 
	      {0,0,0,0}, 
	    };
	
	int refprism_2fa_3e_0v_splitelements[][5] = 
	  {
	      {1,2,3,4,25},
	      {2,1,3,5,26},
	      {3,1,2,6,27}, 
	      {4,1,6,5,37},
	      {5,2,4,6,38},
	      {6,4,5,3,39}, 
	      {0,0,0,0,0},
	  };
	
	  HPREF_ELEMENT_TYPE refprism_2fa_3e_0v_newelstypes[] =
	    {
	      HP_PRISM,
	      HP_HEX,
	      HP_HEX,
	      HP_HEX,
	      HP_PRISM,
	      HP_PRISM,
	      HP_PRISM,
	      HP_PRISM_SINGEDGE,
	      HP_PRISM_SINGEDGE,
	      HP_PRISM_SINGEDGE,
	
	      HP_PRISM_1FA_0E_0V,
	      HP_HEX_1F_0E_0V, 
	      HP_HEX_1F_0E_0V, 
	      HP_HEX_1F_0E_0V, 
	      HP_PRISM_1FA_0E_0V,
	      HP_PRISM_1FA_0E_0V,
	      HP_PRISM_1FA_0E_0V,
	      HP_PRISM_1FA_1E_0V,
	      HP_PRISM_1FA_1E_0V,
	      HP_PRISM_1FA_1E_0V,
	
	      HP_PRISM_1FA_0E_0V,
	      HP_HEX_1F_0E_0V, 
	      HP_HEX_1F_0E_0V, 
	      HP_HEX_1F_0E_0V, 
	      HP_PRISM_1FA_0E_0V,
	      HP_PRISM_1FA_0E_0V,
	      HP_PRISM_1FA_0E_0V,
	      HP_PRISM_1FA_1E_0V,
	      HP_PRISM_1FA_1E_0V,
	      HP_PRISM_1FA_1E_0V,
	
	      HP_NONE
	    };
	
	  int refprism_2fa_3e_0v_newels[][8] =
	    {
	      { 25, 26, 27, 37, 38, 39}, 
	      { 19, 20, 26, 25, 31, 32, 38, 37},  
	      { 27, 26, 21, 22, 39, 38, 33, 34}, 
	      { 23, 24, 25, 27, 35, 36, 37, 39}, 
	      { 19, 25, 24, 31, 37, 36}, 
	      { 26, 20, 21, 38, 32, 33},
	      { 23, 27, 22, 35, 39, 34}, 
	      { 16, 19, 24, 28, 31, 36}, 
	      { 17, 21, 20, 29, 33, 32}, 
	      { 18, 23, 22, 30, 35, 34}, 
	
	      { 13, 14, 15, 25, 26, 27}, 
	      { 7, 8, 14, 13, 19, 20, 26, 25},
	      { 15, 14, 9, 10, 27, 26, 21, 22}, 
	      { 12, 13, 15, 11, 24, 25, 27, 23}, 
	      { 14, 8, 9, 26, 20, 21}, 
	      { 11, 15, 10, 23, 27, 22}, 
	      { 7, 13 , 12, 19, 25, 24}, 
	      { 2, 9, 8, 17, 21, 20}, 
	      { 3, 11, 10, 18, 23, 22}, 
	      { 1, 7, 12, 16, 19, 24}, 
	
	      { 48, 47, 46, 39, 38, 37 }, 
	      { 48, 43, 42, 47, 39, 34, 33, 38}, 
	      { 45, 44, 48, 46, 36, 35, 39, 37},
	      { 46, 47, 41, 40, 37, 38, 32, 31}, 
	      { 47, 42, 41, 38, 33, 32}, 
	      { 45, 46, 40, 36, 37, 31}, 
	      { 44, 43, 48, 35, 34, 39},
	      { 6, 43, 44, 30, 34, 35}, 
	      { 5, 41, 42, 29, 32, 33}, 
	      { 4, 45, 40, 28, 36, 31},
	    };
	
	HPRef_Struct refprism_2fa_3e_0v =
	  {
	    HP_PRISM,
	    refprism_2fa_3e_0v_splitedges, 
	    refprism_2fa_3e_0v_splitfaces, 
	    refprism_2fa_3e_0v_splitelements, 
	    refprism_2fa_3e_0v_newelstypes, 
	    refprism_2fa_3e_0v_newels
	  };
	
	*/
	/*
	
	//  HP_PRISM_1FB_1EB_0V   ... quad face 1-2-4-5
	  int refprism_1fb_1eb_0v_splitedges[][3] =
	    {
	      { 1, 3, 7 },
	      { 2, 3, 8 },
	      { 4, 6, 9 },
	      { 5, 6, 10 },
	      { 2, 1, 11 },
	      { 5, 4, 12 },
	      { 0, 0, 0 }
	    };
	  HPREF_ELEMENT_TYPE refprism_1fb_1eb_0v_newelstypes[] =
	    {
	      HP_HEX_1F_0E_0V,
	      HP_PRISM_1FB_1EB_0V,
	      HP_PRISM,
	      HP_NONE,
	    };
	  int refprism_1fb_1eb_0v_newels[][8] =
	    {
	      { 1, 4, 12, 11, 7, 9, 10, 8  },
	      { 11, 2, 8, 12, 5, 10 },
	      { 7, 8, 3, 9, 10, 6 }
	    };
	  HPRef_Struct refprism_1fb_1eb_0v =
	    {
	      HP_PRISM,
	      refprism_1fb_1eb_0v_splitedges, 
	      0, 0,
	      refprism_1fb_1eb_0v_newelstypes, 
	      refprism_1fb_1eb_0v_newels
	    };
	
	
	
	
	
	
	
	
	
	
	  // HP_PRISM_2F_0E_0V
	  int refprism_2f_0e_0v_splitedges[][3] =
	    {
	      { 1, 3, 7 },
	      { 2, 1, 8 },
	      { 2, 3, 9 },
	      { 3, 1, 10 },
	
	      { 4, 6, 12 },
	      { 5, 4, 13 },
	      { 5, 6, 14 },
	      { 6, 4, 15 },
	
	      { 0, 0, 0 }
	    };
	
	  int refprism_2f_0e_0v_splitfaces[][4] =
	    {
	      { 2, 1, 3, 11 },
	      { 5, 4, 6, 16 },
	      { 0, 0, 0, 0 },
	    };
	
	  HPREF_ELEMENT_TYPE refprism_2f_0e_0v_newelstypes[] =
	    {
	      HP_HEX_1F_0E_0V,
	      HP_HEX_1F_0E_0V,
	      HP_PRISM_1FB_1EA_0V,
	      HP_PRISM_1FB_1EA_0V,
	      HP_PRISM,
	      HP_NONE,
	    };
	  int refprism_2f_0e_0v_newels[][8] =
	    {
	      //{ 1, 8, 11, 7, 4, 13, 16, 12 },
	      // { 9, 3, 10, 11, 14, 6, 15, 16 },
	      { 1, 4, 13, 8, 7, 12, 16, 11 },
	      { 9, 14, 6, 3, 11, 16, 15, 10 },
	      { 2, 9, 11, 5, 14, 16 },
	      // { 8, 2, 11, 13, 5, 16 },
	      { 5, 13, 16, 2, 8, 11 },
	      { 7, 11, 10, 12, 16, 15 }
	    };
	  HPRef_Struct refprism_2f_0e_0v =
	    {
	      HP_PRISM,
	      refprism_2f_0e_0v_splitedges, 
	      refprism_2f_0e_0v_splitfaces, 
	      0,
	      refprism_2f_0e_0v_newelstypes, 
	      refprism_2f_0e_0v_newels
	    };
	
	*/

	  // HP_PYRAMID
	  public static int[][] refpyramid_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refpyramid_newelstypes = {HP_PYRAMID, HP_NONE};
	  public static int[][] refpyramid_newels =
	  {
		  new int[] {1, 2, 3, 4, 5, 0, 0, 0}
	  };
	  public static HPRef_Struct refpyramid = new HPRef_Struct(HP_PYRAMID, refpyramid_splitedges, 0, 0, refpyramid_newelstypes, refpyramid_newels);


	// singular point 1      
	  // HP_PYRAMID_0E_1V
	  public static int[][] refpyramid_0e_1v_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refpyramid_0e_1v_newelstypes = {HP_TET_0E_1V, HP_TET, HP_NONE};
	  public static int[][] refpyramid_0e_1v_newels =
	  {
		  new int[] {1, 2, 4, 5, 0, 0, 0, 0},
		  new int[] {2, 3, 4, 5, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refpyramid_0e_1v = new HPRef_Struct(HP_PYRAMID, refpyramid_0e_1v_splitedges, 0, 0, refpyramid_0e_1v_newelstypes, refpyramid_0e_1v_newels);


	// singular edges 1-2 1-4 singular point 1 
	  // HP_PYRAMID_EDGES
	  public static int[][] refpyramid_edges_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };
	  public static HPREF_ELEMENT_TYPE[] refpyramid_edges_newelstypes = {HP_TET_1E_1VA, HP_TET_1E_1VA, HP_NONE};
	  public static int[][] refpyramid_edges_newels =
	  {
		  new int[] {1, 2, 3, 5, 0, 0, 0, 0},
		  new int[] {1, 4, 5, 3, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refpyramid_edges = new HPRef_Struct(HP_PYRAMID, refpyramid_edges_splitedges, 0, 0, refpyramid_edges_newelstypes, refpyramid_edges_newels);



	// singular face 1-2-5 singular point 5
	  // HP_PYRAMID_1FB_0E_1VA
	  public static int[][] refpyramid_1fb_0e_1va_splitedges =
	  {
		  new int[] {1, 4, 6},
		  new int[] {2, 3, 7},
		  new int[] {5, 1, 8},
		  new int[] {5, 2, 9},
		  new int[] {5, 3, 10},
		  new int[] {5, 4, 11},
		  new int[] {0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refpyramid_1fb_0e_1va_newelstypes = {HP_HEX_1F_0E_0V, HP_PYRAMID_1FB_0E_1VA, HP_PRISM, HP_NONE};
	  public static int[][] refpyramid_1fb_0e_1va_newels =
	  {
		  new int[] {1, 8, 9, 2, 6, 11, 10, 7},
		  new int[] {8, 9, 10, 11, 5, 0, 0, 0},
		  new int[] {3, 7, 10, 4, 6, 11, 0, 0}
	  };
	  public static HPRef_Struct refpyramid_1fb_0e_1va = new HPRef_Struct(HP_PYRAMID, refpyramid_1fb_0e_1va_splitedges, 0, 0, refpyramid_1fb_0e_1va_newelstypes, refpyramid_1fb_0e_1va_newels);





	// HP_QUAD
	public static int[][] refquad_splitedges =
	{
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_newelstypes = {HP_QUAD, HP_NONE};
	public static int[][] refquad_newels =
	{
		new int[] {1, 2, 3, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad = new HPRef_Struct(HP_QUAD, refquad_splitedges, 0, 0, refquad_newelstypes, refquad_newels);







	// HP_QUAD_SINGCORNER
	public static int[][] refquad_singcorner_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_singcorner_newelstypes = {HP_TRIG_SINGCORNER, HP_QUAD, HP_TRIG, HP_NONE};
	public static int[][] refquad_singcorner_newels =
	{
		new int[] {1, 5, 6, 0, 0, 0, 0, 0},
		new int[] {2, 4, 6, 5, 0, 0, 0, 0},
		new int[] {2, 3, 4, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_singcorner = new HPRef_Struct(HP_QUAD, refquad_singcorner_splitedges, 0, 0, refquad_singcorner_newelstypes, refquad_singcorner_newels);





	// HP_DUMMY_QUAD_SINGCORNER
	public static int[][] refdummyquad_singcorner_splitedges =
	{
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refdummyquad_singcorner_newelstypes = {HP_TRIG_SINGCORNER, HP_TRIG, HP_NONE};
	public static int[][] refdummyquad_singcorner_newels =
	{
		new int[] {1, 2, 4, 0, 0, 0, 0, 0},
		new int[] {4, 2, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refdummyquad_singcorner = new HPRef_Struct(HP_QUAD, refdummyquad_singcorner_splitedges, 0, 0, refdummyquad_singcorner_newelstypes, refdummyquad_singcorner_newels);







	// HP_QUAD_SINGEDGE
	public static int[][] refquad_singedge_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_singedge_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_singedge_newels =
	{
		new int[] {1, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 3, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_singedge = new HPRef_Struct(HP_QUAD, refquad_singedge_splitedges, 0, 0, refquad_singedge_newelstypes, refquad_singedge_newels);






	// HP_QUAD_0E_2VA
	public static int[][] refquad_0e_2va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_0e_2va_newelstypes = {HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_0e_2va_newels =
	{
		new int[] {1, 5, 6, 0, 0, 0, 0, 0},
		new int[] {2, 8, 7, 0, 0, 0, 0, 0},
		new int[] {5, 7, 8, 6, 0, 0, 0, 0},
		new int[] {6, 8, 3, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_0e_2va = new HPRef_Struct(HP_QUAD, refquad_0e_2va_splitedges, 0, 0, refquad_0e_2va_newelstypes, refquad_0e_2va_newels);



	// HP_QUAD_0E_2VB
	public static int[][] refquad_0e_2vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {3, 4, 7},
		new int[] {3, 2, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_0e_2vb_newelstypes = {HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_0e_2vb_newels =
	{
		new int[] {1, 5, 6, 0, 0, 0, 0, 0},
		new int[] {3, 7, 8, 0, 0, 0, 0, 0},
		new int[] {5, 2, 4, 6, 0, 0, 0, 0},
		new int[] {2, 8, 7, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_0e_2vb = new HPRef_Struct(HP_QUAD, refquad_0e_2vb_splitedges, 0, 0, refquad_0e_2vb_newelstypes, refquad_0e_2vb_newels);




	// HP_QUAD_0E_3V
	public static int[][] refquad_0e_3v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 2, 9},
		new int[] {3, 4, 10},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_0e_3v_splitfaces =
	{
		new int[] {2, 3, 1, 14},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_0e_3v_newelstypes = {HP_TRIG_SINGCORNER, HP_DUMMY_QUAD_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_0e_3v_newels =
	{
		new int[] {1, 5, 6, 0, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {3, 10, 9, 0, 0, 0, 0, 0},
		new int[] {5, 7, 14, 6, 0, 0, 0, 0},
		new int[] {8, 9, 10, 14, 0, 0, 0, 0},
		new int[] {6, 14, 10, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_0e_3v = new HPRef_Struct(HP_QUAD, refquad_0e_3v_splitedges, refquad_0e_3v_splitfaces, 0, refquad_0e_3v_newelstypes, refquad_0e_3v_newels);




	// HP_QUAD_0E_4V
	public static int[][] refquad_0e_4v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 2, 9},
		new int[] {3, 4, 10},
		new int[] {4, 1, 11},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_0e_4v_splitfaces =
	{
		new int[] {1, 2, 4, 13},
		new int[] {2, 3, 1, 14},
		new int[] {3, 4, 2, 15},
		new int[] {4, 1, 3, 16},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_0e_4v_newelstypes = {HP_DUMMY_QUAD_SINGCORNER, HP_DUMMY_QUAD_SINGCORNER, HP_DUMMY_QUAD_SINGCORNER, HP_DUMMY_QUAD_SINGCORNER, HP_QUAD, HP_QUAD, HP_QUAD, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_0e_4v_newels =
	{
		new int[] {1, 5, 13, 6, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {3, 10, 15, 9, 0, 0, 0, 0},
		new int[] {4, 11, 16, 12, 0, 0, 0, 0},
		new int[] {5, 7, 14, 13, 0, 0, 0, 0},
		new int[] {8, 9, 15, 14, 0, 0, 0, 0},
		new int[] {10, 12, 16, 15, 0, 0, 0, 0},
		new int[] {11, 6, 13, 16, 0, 0, 0, 0},
		new int[] {13, 14, 15, 16, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_0e_4v = new HPRef_Struct(HP_QUAD, refquad_0e_4v_splitedges, refquad_0e_4v_splitfaces, 0, refquad_0e_4v_newelstypes, refquad_0e_4v_newels);








	// HP_QUAD_1E_1VA
	public static int[][] refquad_1e_1va_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_1va_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGEDGECORNER1, HP_NONE};
	public static int[][] refquad_1e_1va_newels =
	{
		new int[] {7, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 3, 4, 0, 0, 0, 0},
		new int[] {1, 7, 5, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_1va = new HPRef_Struct(HP_QUAD, refquad_1e_1va_splitedges, 0, 0, refquad_1e_1va_newelstypes, refquad_1e_1va_newels);




	// HP_QUAD_1E_1VB
	public static int[][] refquad_1e_1vb_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {2, 1, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_1vb_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGEDGECORNER2, HP_NONE};
	public static int[][] refquad_1e_1vb_newels =
	{
		new int[] {1, 7, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 3, 4, 0, 0, 0, 0},
		new int[] {7, 2, 6, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_1vb = new HPRef_Struct(HP_QUAD, refquad_1e_1vb_splitedges, 0, 0, refquad_1e_1vb_newelstypes, refquad_1e_1vb_newels);



	// HP_QUAD_1E_1VC
	public static int[][] refquad_1e_1vc_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {3, 4, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_1vc_newelstypes = {HP_QUAD_SINGEDGE, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_1vc_newels =
	{
		new int[] {1, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 4, 0, 0, 0, 0, 0},
		new int[] {4, 6, 7, 8, 0, 0, 0, 0},
		new int[] {3, 8, 7, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_1vc = new HPRef_Struct(HP_QUAD, refquad_1e_1vc_splitedges, 0, 0, refquad_1e_1vc_newelstypes, refquad_1e_1vc_newels);



	// HP_QUAD_1E_1VD
	public static int[][] refquad_1e_1vd_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {4, 1, 7},
		new int[] {4, 3, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_1vd_newelstypes = {HP_QUAD_SINGEDGE, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_1vd_newels =
	{
		new int[] {1, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 3, 0, 0, 0, 0, 0},
		new int[] {5, 3, 8, 7, 0, 0, 0, 0},
		new int[] {4, 7, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_1vd = new HPRef_Struct(HP_QUAD, refquad_1e_1vd_splitedges, 0, 0, refquad_1e_1vd_newelstypes, refquad_1e_1vd_newels);







	// HP_QUAD_1E_2VA
	public static int[][] refquad_1e_2va_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_2va_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_NONE};
	public static int[][] refquad_1e_2va_newels =
	{
		new int[] {7, 8, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 3, 4, 0, 0, 0, 0},
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {8, 2, 6, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_2va = new HPRef_Struct(HP_QUAD, refquad_1e_2va_splitedges, 0, 0, refquad_1e_2va_newelstypes, refquad_1e_2va_newels);




	// HP_QUAD_1E_2VB
	public static int[][] refquad_1e_2vb_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {3, 2, 8},
		new int[] {3, 4, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_2vb_newelstypes = {HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_2vb_newels =
	{
		new int[] {7, 2, 6, 5, 0, 0, 0, 0},
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {5, 6, 4, 0, 0, 0, 0, 0},
		new int[] {4, 6, 8, 9, 0, 0, 0, 0},
		new int[] {3, 9, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_2vb = new HPRef_Struct(HP_QUAD, refquad_1e_2vb_splitedges, 0, 0, refquad_1e_2vb_newelstypes, refquad_1e_2vb_newels);




	// HP_QUAD_1E_2VC
	public static int[][] refquad_1e_2vc_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {4, 1, 8},
		new int[] {4, 3, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_2vc_newelstypes = {HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_2vc_newels =
	{
		new int[] {7, 2, 6, 5, 0, 0, 0, 0},
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {5, 6, 3, 0, 0, 0, 0, 0},
		new int[] {5, 3, 9, 8, 0, 0, 0, 0},
		new int[] {4, 8, 9, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_2vc = new HPRef_Struct(HP_QUAD, refquad_1e_2vc_splitedges, 0, 0, refquad_1e_2vc_newelstypes, refquad_1e_2vc_newels);




	// HP_QUAD_1E_2VD
	public static int[][] refquad_1e_2vd_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {2, 1, 7},
		new int[] {3, 2, 8},
		new int[] {3, 4, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_2vd_newelstypes = {HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER2, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_2vd_newels =
	{
		new int[] {1, 7, 6, 5, 0, 0, 0, 0},
		new int[] {7, 2, 6, 0, 0, 0, 0, 0},
		new int[] {5, 6, 4, 0, 0, 0, 0, 0},
		new int[] {4, 6, 8, 9, 0, 0, 0, 0},
		new int[] {3, 9, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_2vd = new HPRef_Struct(HP_QUAD, refquad_1e_2vd_splitedges, 0, 0, refquad_1e_2vd_newelstypes, refquad_1e_2vd_newels);





	// HP_QUAD_1E_2VE
	public static int[][] refquad_1e_2ve_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {2, 1, 7},
		new int[] {4, 1, 8},
		new int[] {4, 3, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_2ve_newelstypes = {HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER2, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_2ve_newels =
	{
		new int[] {1, 7, 6, 5, 0, 0, 0, 0},
		new int[] {7, 2, 6, 0, 0, 0, 0, 0},
		new int[] {5, 6, 3, 0, 0, 0, 0, 0},
		new int[] {5, 3, 9, 8, 0, 0, 0, 0},
		new int[] {4, 8, 9, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_2ve = new HPRef_Struct(HP_QUAD, refquad_1e_2ve_splitedges, 0, 0, refquad_1e_2ve_newelstypes, refquad_1e_2ve_newels);






	// HP_QUAD_1E_2VF
	public static int[][] refquad_1e_2vf_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {4, 1, 7},
		new int[] {4, 3, 8},
		new int[] {3, 2, 9},
		new int[] {3, 4, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_2vf_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD, HP_QUAD, HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_1e_2vf_newels =
	{
		new int[] {1, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 9, 7, 0, 0, 0, 0},
		new int[] {7, 9, 10, 8, 0, 0, 0, 0},
		new int[] {4, 7, 8, 0, 0, 0, 0, 0},
		new int[] {3, 10, 9, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_2vf = new HPRef_Struct(HP_QUAD, refquad_1e_2vf_splitedges, 0, 0, refquad_1e_2vf_newelstypes, refquad_1e_2vf_newels);





	// HP_QUAD_1E_3VA
	public static int[][] refquad_1e_3va_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {3, 2, 9},
		new int[] {3, 4, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_3va_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGCORNER, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG, HP_NONE};
	public static int[][] refquad_1e_3va_newels =
	{
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {8, 2, 6, 0, 0, 0, 0, 0},
		new int[] {3, 10, 9, 0, 0, 0, 0, 0},
		new int[] {7, 8, 6, 5, 0, 0, 0, 0},
		new int[] {4, 6, 9, 10, 0, 0, 0, 0},
		new int[] {5, 6, 4, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_3va = new HPRef_Struct(HP_QUAD, refquad_1e_3va_splitedges, 0, 0, refquad_1e_3va_newelstypes, refquad_1e_3va_newels);





	// HP_QUAD_1E_3VB
	public static int[][] refquad_1e_3vb_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {4, 1, 9},
		new int[] {4, 3, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_3vb_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGCORNER, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG, HP_NONE};
	public static int[][] refquad_1e_3vb_newels =
	{
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {8, 2, 6, 0, 0, 0, 0, 0},
		new int[] {4, 9, 10, 0, 0, 0, 0, 0},
		new int[] {7, 8, 6, 5, 0, 0, 0, 0},
		new int[] {5, 3, 10, 9, 0, 0, 0, 0},
		new int[] {5, 6, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_3vb = new HPRef_Struct(HP_QUAD, refquad_1e_3vb_splitedges, 0, 0, refquad_1e_3vb_newelstypes, refquad_1e_3vb_newels);





	// HP_QUAD_1E_3VC
	public static int[][] refquad_1e_3vc_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {3, 2, 8},
		new int[] {3, 4, 9},
		new int[] {4, 3, 10},
		new int[] {4, 1, 11},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_3vc_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD_SINGEDGE, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_1e_3vc_newels =
	{
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {3, 9, 8, 0, 0, 0, 0, 0},
		new int[] {4, 11, 10, 0, 0, 0, 0, 0},
		new int[] {7, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 8, 11, 0, 0, 0, 0},
		new int[] {11, 8, 9, 10, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_3vc = new HPRef_Struct(HP_QUAD, refquad_1e_3vc_splitedges, 0, 0, refquad_1e_3vc_newelstypes, refquad_1e_3vc_newels);




	// HP_QUAD_1E_3VD
	public static int[][] refquad_1e_3vd_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {2, 1, 7},
		new int[] {3, 2, 8},
		new int[] {3, 4, 9},
		new int[] {4, 3, 10},
		new int[] {4, 1, 11},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_3vd_newelstypes = {HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD_SINGEDGE, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_1e_3vd_newels =
	{
		new int[] {7, 2, 6, 0, 0, 0, 0, 0},
		new int[] {3, 9, 8, 0, 0, 0, 0, 0},
		new int[] {4, 11, 10, 0, 0, 0, 0, 0},
		new int[] {1, 7, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 8, 11, 0, 0, 0, 0},
		new int[] {11, 8, 9, 10, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_3vd = new HPRef_Struct(HP_QUAD, refquad_1e_3vd_splitedges, 0, 0, refquad_1e_3vd_newelstypes, refquad_1e_3vd_newels);






	// HP_QUAD_1E_4V
	public static int[][] refquad_1e_4v_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {4, 1, 9},
		new int[] {3, 2, 10},
		new int[] {4, 3, 11},
		new int[] {3, 4, 12},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_1e_4v_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD_SINGEDGE, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] refquad_1e_4v_newels =
	{
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {8, 2, 6, 0, 0, 0, 0, 0},
		new int[] {3, 12, 10, 0, 0, 0, 0, 0},
		new int[] {4, 9, 11, 0, 0, 0, 0, 0},
		new int[] {7, 8, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 10, 9, 0, 0, 0, 0},
		new int[] {9, 10, 12, 11, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_1e_4v = new HPRef_Struct(HP_QUAD, refquad_1e_4v_splitedges, 0, 0, refquad_1e_4v_newelstypes, refquad_1e_4v_newels);

	////////////////////////////////////////////////////////////////////////////////

	// HP_QUAD_2E
	public static int[][] refquad_2e_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 3, 8},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_splitfaces =
	{
		new int[] {1, 2, 4, 9},
		new int[] {0, 0, 0, 0}
	};


	/* 
	   HPREF_ELEMENT_TYPE refquad_2e_newelstypes[] =
	{
	  HP_TRIG_SINGEDGECORNER1,
	  HP_TRIG_SINGEDGECORNER2,
	  HP_QUAD_SINGEDGE,
	  HP_QUAD_SINGEDGE,
	  HP_QUAD,
	  HP_NONE,
	};
	int refquad_2e_newels[][8] =
	{
	  { 1, 5, 9 },
	  { 6, 1, 9 },
	  { 5, 2, 7, 9 },
	  { 4, 6, 9, 8 },
	  { 9, 7, 3, 8 },
	};
	*/ 

	// SZ refine to 4 quads 
	public static HPREF_ELEMENT_TYPE[] refquad_2e_newelstypes = {HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_2e_newels =
	{
		new int[] {1, 5, 9, 6, 0, 0, 0, 0},
		new int[] {5, 2, 7, 9, 0, 0, 0, 0},
		new int[] {4, 6, 9, 8, 0, 0, 0, 0},
		new int[] {9, 7, 3, 8, 0, 0, 0, 0}
	};

	public static HPRef_Struct refquad_2e = new HPRef_Struct(HP_QUAD, refquad_2e_splitedges, refquad_2e_splitfaces, 0, refquad_2e_newelstypes, refquad_2e_newels);


	// HP_QUAD_2E_1VA
	public static int[][] refquad_2e_1va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 3, 8},
		new int[] {2, 1, 10},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_1va_splitfaces =
	{
		new int[] {1, 2, 4, 9},
		new int[] {0, 0, 0, 0}
	};

	/* 
	HPREF_ELEMENT_TYPE refquad_2e_1va_newelstypes[] =
	{
	  HP_TRIG_SINGEDGECORNER1,
	  HP_TRIG_SINGEDGECORNER2,
	  HP_QUAD_SINGEDGE,
	  HP_QUAD_SINGEDGE,
	  HP_QUAD,
	  HP_TRIG_SINGEDGECORNER2,
	  HP_NONE,
	};
	int refquad_2e_1va_newels[][8] =
	{
	  { 1, 5, 9 },
	  { 6, 1, 9 },
	  { 5, 10, 7, 9 },
	  { 4, 6, 9, 8 },
	  { 9, 7, 3, 8 },
	  { 10, 2, 7 },
	};
	*/ 
	// SZ Quad_2e refinement 
	public static HPREF_ELEMENT_TYPE[] refquad_2e_1va_newelstypes = {HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGEDGECORNER2, HP_NONE};
	public static int[][] refquad_2e_1va_newels =
	{
		new int[] {1, 5, 9, 6, 0, 0, 0, 0},
		new int[] {5, 10, 7, 9, 0, 0, 0, 0},
		new int[] {4, 6, 9, 8, 0, 0, 0, 0},
		new int[] {9, 7, 3, 8, 0, 0, 0, 0},
		new int[] {10, 2, 7, 0, 0, 0, 0, 0}
	};

	public static HPRef_Struct refquad_2e_1va = new HPRef_Struct(HP_QUAD, refquad_2e_1va_splitedges, refquad_2e_1va_splitfaces, 0, refquad_2e_1va_newelstypes, refquad_2e_1va_newels);



	// HP_QUAD_2E_1VB
	public static int[][] refquad_2e_1vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 3, 8},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_1vb_splitfaces =
	{
		new int[] {1, 2, 4, 9},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2e_1vb_newelstypes = {HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] refquad_2e_1vb_newels =
	{
		new int[] {1, 5, 9, 6, 0, 0, 0, 0},
		new int[] {5, 2, 7, 9, 0, 0, 0, 0},
		new int[] {4, 6, 9, 8, 0, 0, 0, 0},
		new int[] {7, 8, 9, 0, 0, 0, 0, 0},
		new int[] {8, 7, 10, 11, 0, 0, 0, 0},
		new int[] {3, 11, 10, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2e_1vb = new HPRef_Struct(HP_QUAD, refquad_2e_1vb_splitedges, refquad_2e_1vb_splitfaces, 0, refquad_2e_1vb_newelstypes, refquad_2e_1vb_newels);

	// HP_QUAD_2E_1VC
	public static int[][] refquad_2e_1vc_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 1, 8},
		new int[] {4, 3, 9},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_1vc_splitfaces =
	{
		new int[] {1, 2, 4, 10},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2e_1vc_newelstypes = {HP_QUAD_2E, HP_TRIG_SINGEDGECORNER1, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_2e_1vc_newels =
	{
		new int[] {1, 5, 10, 6, 0, 0, 0, 0},
		new int[] {4, 8, 9, 0, 0, 0, 0, 0},
		new int[] {5, 2, 7, 10, 0, 0, 0, 0},
		new int[] {8, 6, 10, 9, 0, 0, 0, 0},
		new int[] {10, 7, 3, 9, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2e_1vc = new HPRef_Struct(HP_QUAD, refquad_2e_1vc_splitedges, refquad_2e_1vc_splitfaces, 0, refquad_2e_1vc_newelstypes, refquad_2e_1vc_newels);

	// HP_QUAD_2E_2VA
	public static int[][] refquad_2e_2va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 3, 8},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {2, 1, 12},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_2va_splitfaces =
	{
		new int[] {1, 2, 4, 9},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2e_2va_newelstypes = {HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_TRIG_SINGEDGECORNER2, HP_NONE};
	public static int[][] refquad_2e_2va_newels =
	{
		new int[] {1, 5, 9, 6, 0, 0, 0, 0},
		new int[] {5, 12, 7, 9, 0, 0, 0, 0},
		new int[] {4, 6, 9, 8, 0, 0, 0, 0},
		new int[] {7, 8, 9, 0, 0, 0, 0, 0},
		new int[] {8, 7, 10, 11, 0, 0, 0, 0},
		new int[] {3, 11, 10, 0, 0, 0, 0, 0},
		new int[] {12, 2, 7, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2e_2va = new HPRef_Struct(HP_QUAD, refquad_2e_2va_splitedges, refquad_2e_2va_splitfaces, 0, refquad_2e_2va_newelstypes, refquad_2e_2va_newels);






	// HP_QUAD_2E_2VB
	public static int[][] refquad_2e_2vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {4, 1, 9},
		new int[] {4, 3, 10},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_2vb_splitfaces =
	{
		new int[] {1, 2, 4, 11},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2e_2vb_newelstypes = {HP_QUAD_2E, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_2e_2vb_newels =
	{
		new int[] {1, 5, 11, 6, 0, 0, 0, 0},
		new int[] {4, 9, 10, 0, 0, 0, 0, 0},
		new int[] {7, 2, 8, 0, 0, 0, 0, 0},
		new int[] {5, 7, 8, 11, 0, 0, 0, 0},
		new int[] {9, 6, 11, 10, 0, 0, 0, 0},
		new int[] {3, 10, 11, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2e_2vb = new HPRef_Struct(HP_QUAD, refquad_2e_2vb_splitedges, refquad_2e_2vb_splitfaces, 0, refquad_2e_2vb_newelstypes, refquad_2e_2vb_newels);

	// HP_QUAD_2E_2VC
	public static int[][] refquad_2e_2vc_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 3, 8},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {4, 1, 12},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_2vc_splitfaces =
	{
		new int[] {1, 2, 4, 9},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2e_2vc_newelstypes = {HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_TRIG_SINGEDGECORNER1, HP_NONE};
	public static int[][] refquad_2e_2vc_newels =
	{
		new int[] {1, 5, 9, 0, 0, 0, 0, 0},
		new int[] {6, 1, 9, 0, 0, 0, 0, 0},
		new int[] {5, 2, 7, 9, 0, 0, 0, 0},
		new int[] {12, 6, 9, 8, 0, 0, 0, 0},
		new int[] {7, 8, 9, 0, 0, 0, 0, 0},
		new int[] {8, 7, 10, 11, 0, 0, 0, 0},
		new int[] {3, 11, 10, 0, 0, 0, 0, 0},
		new int[] {4, 12, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2e_2vc = new HPRef_Struct(HP_QUAD, refquad_2e_2vc_splitedges, refquad_2e_2vc_splitfaces, 0, refquad_2e_2vc_newelstypes, refquad_2e_2vc_newels);

	// HP_QUAD_2E_3V  
	public static int[][] refquad_2e_3v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {4, 3, 8},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {2, 1, 12},
		new int[] {4, 1, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2e_3v_splitfaces =
	{
		new int[] {1, 2, 4, 9},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2e_3v_newelstypes = {HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG, HP_QUAD, HP_TRIG_SINGCORNER, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_NONE};
	public static int[][] refquad_2e_3v_newels =
	{
		new int[] {1, 5, 9, 6, 0, 0, 0, 0},
		new int[] {5, 12, 7, 9, 0, 0, 0, 0},
		new int[] {13, 6, 9, 8, 0, 0, 0, 0},
		new int[] {7, 8, 9, 0, 0, 0, 0, 0},
		new int[] {8, 7, 10, 11, 0, 0, 0, 0},
		new int[] {3, 11, 10, 0, 0, 0, 0, 0},
		new int[] {12, 2, 7, 0, 0, 0, 0, 0},
		new int[] {4, 13, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2e_3v = new HPRef_Struct(HP_QUAD, refquad_2e_3v_splitedges, refquad_2e_3v_splitfaces, 0, refquad_2e_3v_newelstypes, refquad_2e_3v_newels);

	// HP_QUAD_2EB_0V
	public static int[][] refquad_2eb_0v_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2eb_0v_splitfaces =
	{
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_0v_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_0v_newels =
	{
		new int[] {1, 2, 6, 5, 0, 0, 0, 0},
		new int[] {3, 4, 8, 7, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_0v = new HPRef_Struct(HP_QUAD, refquad_2eb_0v_splitedges, refquad_2eb_0v_splitfaces, 0, refquad_2eb_0v_newelstypes, refquad_2eb_0v_newels);


	// HP_QUAD_2EB_1VA
	public static int[][] refquad_2eb_1va_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {1, 2, 9},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2eb_1va_splitfaces =
	{
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_1va_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_1va_newels =
	{
		new int[] {9, 2, 6, 5, 0, 0, 0, 0},
		new int[] {3, 4, 8, 7, 0, 0, 0, 0},
		new int[] {1, 9, 5, 0, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_1va = new HPRef_Struct(HP_QUAD, refquad_2eb_1va_splitedges, refquad_2eb_1va_splitfaces, 0, refquad_2eb_1va_newelstypes, refquad_2eb_1va_newels);

	// HP_QUAD_2EB_1VB
	public static int[][] refquad_2eb_1vb_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {2, 1, 9},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2eb_1vb_splitfaces =
	{
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_1vb_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER2, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_1vb_newels =
	{
		new int[] {1, 9, 6, 5, 0, 0, 0, 0},
		new int[] {3, 4, 8, 7, 0, 0, 0, 0},
		new int[] {9, 2, 6, 0, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_1vb = new HPRef_Struct(HP_QUAD, refquad_2eb_1vb_splitedges, refquad_2eb_1vb_splitfaces, 0, refquad_2eb_1vb_newelstypes, refquad_2eb_1vb_newels);

	// HP_QUAD_2EB_2VA
	public static int[][] refquad_2eb_2va_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {1, 2, 9},
		new int[] {2, 1, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_2va_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_2va_newels =
	{
		new int[] {9, 10, 6, 5, 0, 0, 0, 0},
		new int[] {3, 4, 8, 7, 0, 0, 0, 0},
		new int[] {1, 9, 5, 0, 0, 0, 0, 0},
		new int[] {10, 2, 6, 0, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_2va = new HPRef_Struct(HP_QUAD, refquad_2eb_2va_splitedges, 0, 0, refquad_2eb_2va_newelstypes, refquad_2eb_2va_newels);



	// HP_QUAD_2EB_2VB
	public static int[][] refquad_2eb_2vb_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {1, 2, 9},
		new int[] {3, 4, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_2vb_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER1, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_2vb_newels =
	{
		new int[] {9, 2, 6, 5, 0, 0, 0, 0},
		new int[] {10, 4, 8, 7, 0, 0, 0, 0},
		new int[] {1, 9, 5, 0, 0, 0, 0, 0},
		new int[] {3, 10, 7, 0, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_2vb = new HPRef_Struct(HP_QUAD, refquad_2eb_2vb_splitedges, 0, 0, refquad_2eb_2vb_newelstypes, refquad_2eb_2vb_newels);



	// HP_QUAD_2EB_2VC
	public static int[][] refquad_2eb_2vc_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {1, 2, 9},
		new int[] {4, 3, 10},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2eb_2vc_splitfaces =
	{
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_2vc_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_2vc_newels =
	{
		new int[] {9, 2, 6, 5, 0, 0, 0, 0},
		new int[] {3, 10, 8, 7, 0, 0, 0, 0},
		new int[] {1, 9, 5, 0, 0, 0, 0, 0},
		new int[] {10, 4, 8, 0, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_2vc = new HPRef_Struct(HP_QUAD, refquad_2eb_2vc_splitedges, refquad_2eb_2vc_splitfaces, 0, refquad_2eb_2vc_newelstypes, refquad_2eb_2vc_newels);


	// HP_QUAD_2EB_2VD
	public static int[][] refquad_2eb_2vd_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {2, 1, 9},
		new int[] {4, 3, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_2vd_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER2, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_2vd_newels =
	{
		new int[] {1, 9, 6, 5, 0, 0, 0, 0},
		new int[] {3, 10, 8, 7, 0, 0, 0, 0},
		new int[] {9, 2, 6, 0, 0, 0, 0, 0},
		new int[] {10, 4, 8, 0, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_2vd = new HPRef_Struct(HP_QUAD, refquad_2eb_2vd_splitedges, 0, 0, refquad_2eb_2vd_newelstypes, refquad_2eb_2vd_newels);


	// HP_QUAD_2EB_3VA
	public static int[][] refquad_2eb_3va_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {3, 2, 9},
		new int[] {4, 1, 10},
		new int[] {3, 4, 11},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_3va_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_3va_newels =
	{
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {8, 2, 6, 0, 0, 0, 0, 0},
		new int[] {3, 11, 9, 0, 0, 0, 0, 0},
		new int[] {7, 8, 6, 5, 0, 0, 0, 0},
		new int[] {11, 4, 10, 9, 0, 0, 0, 0},
		new int[] {5, 6, 9, 10, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_3va = new HPRef_Struct(HP_QUAD, refquad_2eb_3va_splitedges, 0, 0, refquad_2eb_3va_newelstypes, refquad_2eb_3va_newels);


	// HP_QUAD_2EB_3VB
	public static int[][] refquad_2eb_3vb_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {1, 2, 7},
		new int[] {2, 1, 8},
		new int[] {3, 2, 9},
		new int[] {4, 1, 10},
		new int[] {4, 3, 11},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_3vb_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_2eb_3vb_newels =
	{
		new int[] {1, 7, 5, 0, 0, 0, 0, 0},
		new int[] {8, 2, 6, 0, 0, 0, 0, 0},
		new int[] {11, 4, 10, 0, 0, 0, 0, 0},
		new int[] {7, 8, 6, 5, 0, 0, 0, 0},
		new int[] {3, 11, 10, 9, 0, 0, 0, 0},
		new int[] {5, 6, 9, 10, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_3vb = new HPRef_Struct(HP_QUAD, refquad_2eb_3vb_splitedges, 0, 0, refquad_2eb_3vb_newelstypes, refquad_2eb_3vb_newels);


	// HP_QUAD_2EB_4V
	public static int[][] refquad_2eb_4v_splitedges =
	{
		new int[] {1, 4, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {4, 1, 8},
		new int[] {1, 2, 9},
		new int[] {2, 1, 10},
		new int[] {3, 4, 11},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};
	public static int[][] refquad_2eb_4v_splitfaces =
	{
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] refquad_2eb_4v_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_NONE};
	public static int[][] refquad_2eb_4v_newels =
	{
		new int[] {9, 10, 6, 5, 0, 0, 0, 0},
		new int[] {11, 12, 8, 7, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 0, 0, 0, 0},
		new int[] {1, 9, 5, 0, 0, 0, 0, 0},
		new int[] {10, 2, 6, 0, 0, 0, 0, 0},
		new int[] {3, 11, 7, 0, 0, 0, 0, 0},
		new int[] {12, 4, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_2eb_4v = new HPRef_Struct(HP_QUAD, refquad_2eb_4v_splitedges, refquad_2eb_4v_splitfaces, 0, refquad_2eb_4v_newelstypes, refquad_2eb_4v_newels);



	// HP_QUAD_3E
	public static int[][] refquad_3e_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 4, 10},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_3e_splitfaces =
	{
		new int[] {1, 2, 4, 13},
		new int[] {2, 3, 1, 14},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_3e_newelstypes = {HP_QUAD_2E, HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_3e_newels =
	{
		new int[] {1, 5, 13, 6, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {5, 7, 14, 13, 0, 0, 0, 0},
		new int[] {8, 3, 10, 14, 0, 0, 0, 0},
		new int[] {4, 6, 13, 12, 0, 0, 0, 0},
		new int[] {13, 14, 10, 12, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_3e = new HPRef_Struct(HP_QUAD, refquad_3e_splitedges, refquad_3e_splitfaces, 0, refquad_3e_newelstypes, refquad_3e_newels);







	// HP_QUAD_3E_3VA
	public static int[][] refquad_3e_3va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 4, 10},
		new int[] {3, 2, 11},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_3e_3va_splitfaces =
	{
		new int[] {1, 2, 4, 13},
		new int[] {2, 3, 1, 14},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_3e_3va_newelstypes = {HP_QUAD_2E, HP_QUAD_2E, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_3e_3va_newels =
	{
		new int[] {1, 5, 13, 6, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {11, 3, 10, 0, 0, 0, 0, 0},
		new int[] {5, 7, 14, 13, 0, 0, 0, 0},
		new int[] {8, 11, 10, 14, 0, 0, 0, 0},
		new int[] {4, 6, 13, 12, 0, 0, 0, 0},
		new int[] {13, 14, 10, 12, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_3e_3va = new HPRef_Struct(HP_QUAD, refquad_3e_3va_splitedges, refquad_3e_3va_splitfaces, 0, refquad_3e_3va_newelstypes, refquad_3e_3va_newels);

	// HP_QUAD_3E_3VB
	public static int[][] refquad_3e_3vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 4, 10},
		new int[] {4, 1, 11},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_3e_3vb_splitfaces =
	{
		new int[] {1, 2, 4, 13},
		new int[] {2, 3, 1, 14},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_3e_3vb_newelstypes = {HP_QUAD_2E, HP_QUAD_2E, HP_TRIG_SINGEDGECORNER1, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_3e_3vb_newels =
	{
		new int[] {1, 5, 13, 6, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {4, 11, 12, 0, 0, 0, 0, 0},
		new int[] {5, 7, 14, 13, 0, 0, 0, 0},
		new int[] {8, 3, 10, 14, 0, 0, 0, 0},
		new int[] {11, 6, 13, 12, 0, 0, 0, 0},
		new int[] {13, 14, 10, 12, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_3e_3vb = new HPRef_Struct(HP_QUAD, refquad_3e_3vb_splitedges, refquad_3e_3vb_splitfaces, 0, refquad_3e_3vb_newelstypes, refquad_3e_3vb_newels);









	// HP_QUAD_3E_4V
	public static int[][] refquad_3e_4v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 4, 10},
		new int[] {3, 2, 11},
		new int[] {4, 3, 12},
		new int[] {4, 1, 15},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_3e_4v_splitfaces =
	{
		new int[] {1, 2, 4, 13},
		new int[] {2, 3, 1, 14},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_3e_4v_newelstypes = {HP_QUAD_2E, HP_QUAD_2E, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_3e_4v_newels =
	{
		new int[] {1, 5, 13, 6, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {11, 3, 10, 0, 0, 0, 0, 0},
		new int[] {4, 15, 12, 0, 0, 0, 0, 0},
		new int[] {5, 7, 14, 13, 0, 0, 0, 0},
		new int[] {8, 11, 10, 14, 0, 0, 0, 0},
		new int[] {15, 6, 13, 12, 0, 0, 0, 0},
		new int[] {13, 14, 10, 12, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_3e_4v = new HPRef_Struct(HP_QUAD, refquad_3e_4v_splitedges, refquad_3e_4v_splitfaces, 0, refquad_3e_4v_newelstypes, refquad_3e_4v_newels);









	// HP_QUAD_4E
	public static int[][] refquad_4e_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {3, 2, 9},
		new int[] {3, 4, 10},
		new int[] {4, 1, 11},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};

	public static int[][] refquad_4e_splitfaces =
	{
		new int[] {1, 2, 4, 13},
		new int[] {2, 3, 1, 14},
		new int[] {3, 4, 2, 15},
		new int[] {4, 1, 3, 16},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] refquad_4e_newelstypes = {HP_QUAD_2E, HP_QUAD_2E, HP_QUAD_2E, HP_QUAD_2E, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD, HP_NONE};
	public static int[][] refquad_4e_newels =
	{
		new int[] {1, 5, 13, 6, 0, 0, 0, 0},
		new int[] {2, 8, 14, 7, 0, 0, 0, 0},
		new int[] {3, 10, 15, 9, 0, 0, 0, 0},
		new int[] {4, 11, 16, 12, 0, 0, 0, 0},
		new int[] {5, 7, 14, 13, 0, 0, 0, 0},
		new int[] {8, 9, 15, 14, 0, 0, 0, 0},
		new int[] {10, 12, 16, 15, 0, 0, 0, 0},
		new int[] {11, 6, 13, 16, 0, 0, 0, 0},
		new int[] {13, 14, 15, 16, 0, 0, 0, 0}
	};
	public static HPRef_Struct refquad_4e = new HPRef_Struct(HP_QUAD, refquad_4e_splitedges, refquad_4e_splitfaces, 0, refquad_4e_newelstypes, refquad_4e_newels);

	  // HP_SEGM
	  public static int[][] refsegm_splitedges =
	  {
		  new int[] {0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refsegm_newelstypes = {HP_SEGM, HP_NONE};
	  public static int[][] refsegm_newels =
	  {
		  new int[] {1, 2, 0, 0, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refsegm = new HPRef_Struct(HP_SEGM, refsegm_splitedges, 0, 0, refsegm_newelstypes, refsegm_newels);

	  // HP_SEGM_SINGCORNERL = 2,
	  public static int[][] refsegm_scl_splitedges =
	  {
		  new int[] {1, 2, 3},
		  new int[] {0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refsegm_scl_newelstypes = {HP_SEGM_SINGCORNERL, HP_SEGM, HP_NONE};

	  public static int[][] refsegm_scl_newels =
	  {
		  new int[] {1, 3, 0, 0, 0, 0, 0, 0},
		  new int[] {3, 2, 0, 0, 0, 0, 0, 0},
		  new int[] {0, 0, 0, 0, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refsegm_scl = new HPRef_Struct(HP_SEGM, refsegm_scl_splitedges, 0, 0, refsegm_scl_newelstypes, refsegm_scl_newels);



	  // HP_SEGM_SINGCORNERR
	  public static int[][] refsegm_scr_splitedges =
	  {
		  new int[] {2, 1, 3},
		  new int[] {0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refsegm_scr_newelstypes = {HP_SEGM, HP_SEGM_SINGCORNERR, HP_NONE};
	  public static int[][] refsegm_scr_newels =
	  {
		  new int[] {1, 3, 0, 0, 0, 0, 0, 0},
		  new int[] {3, 2, 0, 0, 0, 0, 0, 0},
		  new int[] {0, 0, 0, 0, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refsegm_scr = new HPRef_Struct(HP_SEGM, refsegm_scr_splitedges, 0, 0, refsegm_scr_newelstypes, refsegm_scr_newels);






	  // HP_SEGM_SINGCORNERS = 3,
	  public static int[][] refsegm_sc2_splitedges =
	  {
		  new int[] {1, 2, 3},
		  new int[] {2, 1, 4},
		  new int[] {0, 0, 0}
	  };

	  public static HPREF_ELEMENT_TYPE[] refsegm_sc2_newelstypes = {HP_SEGM_SINGCORNERL, HP_SEGM_SINGCORNERR, HP_SEGM, HP_NONE};
	  public static int[][] refsegm_sc2_newels =
	  {
		  new int[] {1, 3, 0, 0, 0, 0, 0, 0},
		  new int[] {4, 2, 0, 0, 0, 0, 0, 0},
		  new int[] {3, 4, 0, 0, 0, 0, 0, 0},
		  new int[] {0, 0, 0, 0, 0, 0, 0, 0}
	  };
	  public static HPRef_Struct refsegm_sc2 = new HPRef_Struct(HP_SEGM, refsegm_sc2_splitedges, 0, 0, refsegm_sc2_newelstypes, refsegm_sc2_newels);





	// HP_TET
	public static int[][] reftet_splitedges =
	{
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_newelstypes = {HP_TET, HP_NONE};
	public static int[][] reftet_newels =
	{
		new int[] {1, 2, 3, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet = new HPRef_Struct(HP_TET, reftet_splitedges, 0, 0, reftet_newelstypes, reftet_newels);



	/* *********** Tet - Refinement - 0 edges *************** */

	// HP_TET_0E_1V
	public static int[][] reftet_0e_1v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_0e_1v_newelstypes = {HP_TET_0E_1V, HP_PRISM, HP_NONE};
	public static int[][] reftet_0e_1v_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {5, 6, 7, 2, 3, 4, 0, 0}
	};
	public static HPRef_Struct reftet_0e_1v = new HPRef_Struct(HP_TET, reftet_0e_1v_splitedges, 0, 0, reftet_0e_1v_newelstypes, reftet_0e_1v_newels);



	// HP_TET_0E_2V
	public static int[][] reftet_0e_2v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_0e_2v_newelstypes = {HP_TET_0E_1V, HP_TET_0E_1V, HP_PRISM, HP_PRISM, HP_NONE};
	public static int[][] reftet_0e_2v_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 10, 9, 8, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {4, 10, 7, 3, 9, 6, 0, 0}
	};
	public static HPRef_Struct reftet_0e_2v = new HPRef_Struct(HP_TET, reftet_0e_2v_splitedges, 0, 0, reftet_0e_2v_newelstypes, reftet_0e_2v_newels);





	// HP_TET_0E_3V
	public static int[][] reftet_0e_3v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_0e_3v_splitfaces =
	{
		new int[] {1, 2, 3, 14},
		new int[] {2, 3, 1, 15},
		new int[] {3, 1, 2, 16},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_0e_3v_newelstypes = {HP_PYRAMID_0E_1V, HP_PYRAMID_0E_1V, HP_PYRAMID_0E_1V, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_0e_3v_newels =
	{
		new int[] {1, 5, 14, 6, 7, 0, 0, 0},
		new int[] {2, 9, 15, 8, 10, 0, 0, 0},
		new int[] {3, 11, 16, 12, 13, 0, 0, 0},
		new int[] {5, 14, 7, 8, 15, 10, 0, 0},
		new int[] {9, 15, 10, 12, 16, 13, 0, 0},
		new int[] {6, 7, 14, 11, 13, 16, 0, 0},
		new int[] {14, 15, 16, 7, 10, 13, 0, 0},
		new int[] {7, 10, 13, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_0e_3v = new HPRef_Struct(HP_TET, reftet_0e_3v_splitedges, reftet_0e_3v_splitfaces, 0, reftet_0e_3v_newelstypes, reftet_0e_3v_newels);





	// HP_TET_0E_4V
	public static int[][] reftet_0e_4v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_0e_4v_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {1, 2, 4, 18},
		new int[] {1, 3, 4, 19},
		new int[] {2, 1, 3, 20},
		new int[] {2, 1, 4, 21},
		new int[] {2, 3, 4, 22},
		new int[] {3, 1, 2, 23},
		new int[] {3, 1, 4, 24},
		new int[] {3, 2, 4, 25},
		new int[] {4, 1, 2, 26},
		new int[] {4, 1, 3, 27},
		new int[] {4, 2, 3, 28},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_0e_4v_splitelements =
	{
		new int[] {1, 2, 3, 4, 29},
		new int[] {2, 3, 4, 1, 30},
		new int[] {3, 4, 1, 2, 31},
		new int[] {4, 1, 2, 3, 32},
		new int[] {0, 0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_0e_4v_newelstypes = {HP_HEX_0E_1V, HP_HEX_0E_1V, HP_HEX_0E_1V, HP_HEX_0E_1V, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_0e_4v_newels =
	{
		new int[] {1, 5, 17, 6, 7, 18, 29, 19},
		new int[] {2, 9, 20, 8, 10, 22, 30, 21},
		new int[] {3, 11, 23, 12, 13, 24, 31, 25},
		new int[] {4, 15, 26, 14, 16, 28, 32, 27},
		new int[] {5, 17, 18, 8, 20, 21, 0, 0},
		new int[] {18, 17, 29, 21, 20, 30, 0, 0},
		new int[] {6, 19, 17, 11, 24, 23, 0, 0},
		new int[] {17, 19, 29, 23, 24, 31, 0, 0},
		new int[] {7, 18, 19, 14, 26, 27, 0, 0},
		new int[] {19, 18, 29, 27, 26, 32, 0, 0},
		new int[] {9, 20, 22, 12, 23, 25, 0, 0},
		new int[] {22, 20, 30, 25, 23, 31, 0, 0},
		new int[] {10, 22, 21, 15, 28, 26, 0, 0},
		new int[] {21, 22, 30, 26, 28, 32, 0, 0},
		new int[] {13, 24, 25, 16, 27, 28, 0, 0},
		new int[] {25, 24, 31, 28, 27, 32, 0, 0},
		new int[] {17, 20, 23, 29, 30, 31, 0, 0},
		new int[] {18, 26, 21, 29, 32, 30, 0, 0},
		new int[] {19, 24, 27, 29, 31, 32, 0, 0},
		new int[] {22, 28, 25, 30, 32, 31, 0, 0},
		new int[] {29, 30, 31, 32, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_0e_4v = new HPRef_Struct(HP_TET, reftet_0e_4v_splitedges, reftet_0e_4v_splitfaces, reftet_0e_4v_splitelements, reftet_0e_4v_newelstypes, reftet_0e_4v_newels);

















	/* *********** Tet - Refinement - 1 edge *************** */



	// HP_TET_1E_0V
	public static int[][] reftet_1e_0v_splitedges =
	{
		new int[] {1, 3, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {2, 4, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1e_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM, HP_NONE};
	public static int[][] reftet_1e_0v_newels =
	{
		new int[] {1, 5, 6, 2, 7, 8, 0, 0},
		new int[] {7, 3, 5, 8, 4, 6, 0, 0}
	};
	public static HPRef_Struct reftet_1e_0v = new HPRef_Struct(HP_TET, reftet_1e_0v_splitedges, 0, 0, reftet_1e_0v_newelstypes, reftet_1e_0v_newels);





	// HP_TET_1E_1VA
	public static int[][] reftet_1e_1va_splitedges =
	{
		new int[] {1, 3, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {2, 4, 8},
		new int[] {1, 2, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1e_1va_newelstypes = {HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_NONE};
	public static int[][] reftet_1e_1va_newels =
	{
		new int[] {1, 9, 5, 6, 0, 0, 0, 0},
		new int[] {9, 5, 6, 2, 7, 8, 0, 0},
		new int[] {7, 3, 5, 8, 4, 6, 0, 0}
	};
	public static HPRef_Struct reftet_1e_1va = new HPRef_Struct(HP_TET, reftet_1e_1va_splitedges, 0, 0, reftet_1e_1va_newelstypes, reftet_1e_1va_newels);






	// HP_TET_1E_1VB
	public static int[][] reftet_1e_1vb_splitedges =
	{
		new int[] {1, 3, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {2, 4, 8},
		new int[] {4, 1, 9},
		new int[] {4, 2, 10},
		new int[] {4, 3, 11},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_1vb_splitelements =
	{
		new int[] {4, 1, 2, 3, 12},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_1e_1vb_newelstypes = {HP_PRISM_SINGEDGE, HP_TET_0E_1V, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_1e_1vb_newels =
	{
		new int[] {1, 5, 6, 2, 7, 8, 0, 0},
		new int[] {4, 11, 10, 9, 0, 0, 0, 0},
		new int[] {7, 8, 10, 11, 12, 0, 0, 0},
		new int[] {3, 7, 11, 12, 0, 0, 0, 0},
		new int[] {5, 11, 9, 6, 12, 0, 0, 0},
		new int[] {5, 3, 11, 12, 0, 0, 0, 0},
		new int[] {6, 9, 10, 8, 12, 0, 0, 0},
		new int[] {5, 7, 3, 12, 0, 0, 0, 0},
		new int[] {5, 6, 8, 7, 12, 0, 0, 0},
		new int[] {9, 11, 10, 12, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1e_1vb = new HPRef_Struct(HP_TET, reftet_1e_1vb_splitedges, 0, reftet_1e_1vb_splitelements, reftet_1e_1vb_newelstypes, reftet_1e_1vb_newels);








	// HP_TET_1E_2VA
	public static int[][] reftet_1e_2va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1e_2va_newelstypes = {HP_TET_1E_1VA, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_NONE};
	public static int[][] reftet_1e_2va_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {4, 10, 7, 3, 9, 6, 0, 0}
	};
	public static HPRef_Struct reftet_1e_2va = new HPRef_Struct(HP_TET, reftet_1e_2va_splitedges, 0, 0, reftet_1e_2va_newelstypes, reftet_1e_2va_newels);







	// HP_TET_1E_2VB
	public static int[][] reftet_1e_2vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {3, 1, 10},
		new int[] {3, 2, 11},
		new int[] {3, 4, 12},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_2vb_splitelements =
	{
		new int[] {3, 4, 1, 2, 13},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_1e_2vb_newelstypes = {HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_TET_0E_1V, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_1e_2vb_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {5, 6, 7, 2, 8, 9, 0, 0},
		new int[] {3, 10, 11, 12, 0, 0, 0, 0},
		new int[] {8, 9, 12, 11, 13, 0, 0, 0},
		new int[] {4, 12, 9, 13, 0, 0, 0, 0},
		new int[] {6, 10, 12, 7, 13, 0, 0, 0},
		new int[] {4, 7, 12, 13, 0, 0, 0, 0},
		new int[] {6, 8, 11, 10, 13, 0, 0, 0},
		new int[] {4, 9, 7, 13, 0, 0, 0, 0},
		new int[] {6, 7, 9, 8, 13, 0, 0, 0},
		new int[] {10, 11, 12, 13, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1e_2vb = new HPRef_Struct(HP_TET, reftet_1e_2vb_splitedges, 0, reftet_1e_2vb_splitelements, reftet_1e_2vb_newelstypes, reftet_1e_2vb_newels);






	// HP_TET_1E_2VC
	public static int[][] reftet_1e_2vc_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {4, 1, 10},
		new int[] {4, 2, 11},
		new int[] {4, 3, 12},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_2vc_splitelements =
	{
		new int[] {4, 1, 2, 3, 13},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_1e_2vc_newelstypes = {HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_TET_0E_1V, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_1e_2vc_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {5, 6, 7, 2, 8, 9, 0, 0},
		new int[] {4, 11, 10, 12, 0, 0, 0, 0},
		new int[] {8, 9, 11, 12, 13, 0, 0, 0},
		new int[] {3, 8, 12, 13, 0, 0, 0, 0},
		new int[] {7, 6, 12, 10, 13, 0, 0, 0},
		new int[] {3, 12, 6, 13, 0, 0, 0, 0},
		new int[] {9, 7, 10, 11, 13, 0, 0, 0},
		new int[] {3, 6, 8, 13, 0, 0, 0, 0},
		new int[] {6, 7, 9, 8, 13, 0, 0, 0},
		new int[] {10, 12, 11, 13, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1e_2vc = new HPRef_Struct(HP_TET, reftet_1e_2vc_splitedges, 0, reftet_1e_2vc_splitelements, reftet_1e_2vc_newelstypes, reftet_1e_2vc_newels);








	/*
	
	// HP_TET_1E_2VD
	int reftet_1e_2vd_splitedges[][3] =
	{
	  { 1, 3, 5 },
	  { 1, 4, 6 },
	  { 2, 3, 7 },
	  { 2, 4, 8 },
	  { 3, 1, 9 },
	  { 3, 2, 10 },
	  { 3, 4, 11 },
	  { 4, 1, 12 },
	  { 4, 2, 13 },
	  { 4, 3, 14 },
	  { 0, 0, 0 }
	};
	HPREF_ELEMENT_TYPE reftet_1e_2vd_newelstypes[] =
	{
	  HP_PRISM_SINGEDGE,
	  HP_TET_0E_1V,
	  HP_TET_0E_1V,
	  HP_PRISM,
	  HP_HEX,
	  HP_NONE,
	};
	int reftet_1e_2vd_newels[][8] =
	{
	  { 1, 5, 6, 2, 7, 8 },
	  { 4, 13, 12, 14 },
	  { 3, 10, 11, 9 },
	  { 14, 13, 12, 11, 10, 9 },
	  { 6, 12, 13, 8, 5, 9, 10, 7 },
	};
	HPRef_Struct reftet_1e_2vd =
	{
	  HP_TET,
	  reftet_1e_2vd_splitedges, 
	  0, 0,
	  reftet_1e_2vd_newelstypes, 
	  reftet_1e_2vd_newels
	};
	
	*/




	//  HP_TET_1E_2VD,  // 1 v on edge
	public static int[][] reftet_1e_2vd_splitedges =
	{
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_2vd_splitfaces =
	{
		new int[] {1, 3, 4, 19},
		new int[] {2, 3, 4, 22},
		new int[] {3, 1, 4, 24},
		new int[] {3, 2, 4, 25},
		new int[] {4, 1, 3, 27},
		new int[] {4, 2, 3, 28},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1e_2vd_newelstypes = {HP_PRISM_SINGEDGE, HP_TET_0E_1V, HP_TET_0E_1V, HP_PRISM, HP_HEX, HP_PYRAMID, HP_HEX, HP_PYRAMID, HP_PRISM, HP_PRISM, HP_NONE};
	public static int[][] reftet_1e_2vd_newels =
	{
		new int[] {1, 6, 7, 2, 9, 10, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {4, 16, 15, 14, 0, 0, 0, 0},
		new int[] {7, 6, 19, 10, 9, 22, 0, 0},
		new int[] {7, 19, 27, 14, 10, 22, 28, 15},
		new int[] {14, 15, 28, 27, 16, 0, 0, 0},
		new int[] {9, 6, 19, 22, 12, 11, 24, 25},
		new int[] {12, 11, 24, 25, 13, 0, 0, 0},
		new int[] {19, 24, 27, 22, 25, 28, 0, 0},
		new int[] {16, 28, 27, 13, 25, 24, 0, 0}
	};
	public static HPRef_Struct reftet_1e_2vd = new HPRef_Struct(HP_TET, reftet_1e_2vd_splitedges, reftet_1e_2vd_splitfaces, 0, reftet_1e_2vd_newelstypes, reftet_1e_2vd_newels);















	// HP_TET_1E_3VA
	public static int[][] reftet_1e_3va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_3va_splitelements =
	{
		new int[] {1, 2, 3, 4, 14},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_1e_3va_newelstypes = {HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_0E_1V, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_1e_3va_newels =
	{
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {6, 7, 10, 9, 14, 0, 0, 0},
		new int[] {4, 10, 7, 14, 0, 0, 0, 0},
		new int[] {9, 10, 13, 12, 14, 0, 0, 0},
		new int[] {4, 13, 10, 14, 0, 0, 0, 0},
		new int[] {6, 11, 13, 7, 14, 0, 0, 0},
		new int[] {4, 7, 13, 14, 0, 0, 0, 0},
		new int[] {6, 11, 12, 9, 14, 0, 0, 0},
		new int[] {11, 13, 12, 14, 0, 0, 0, 0}
	};

	public static HPRef_Struct reftet_1e_3va = new HPRef_Struct(HP_TET, reftet_1e_3va_splitedges, 0, reftet_1e_3va_splitelements, reftet_1e_3va_newelstypes, reftet_1e_3va_newels);






















	//  HP_TET_1E_3VB,  // 1 v on edge
	public static int[][] reftet_1e_3vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_3vb_splitfaces =
	{
		new int[] {1, 3, 4, 19},
		new int[] {2, 3, 4, 22},
		new int[] {3, 1, 4, 24},
		new int[] {3, 2, 4, 25},
		new int[] {4, 1, 3, 27},
		new int[] {4, 2, 3, 28},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1e_3vb_newelstypes = {HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_TET_0E_1V, HP_TET_0E_1V, HP_PRISM, HP_HEX, HP_PYRAMID, HP_HEX, HP_PYRAMID, HP_PRISM, HP_PRISM, HP_NONE};
	public static int[][] reftet_1e_3vb_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {5, 6, 7, 2, 9, 10, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {4, 16, 15, 14, 0, 0, 0, 0},
		new int[] {7, 6, 19, 10, 9, 22, 0, 0},
		new int[] {7, 19, 27, 14, 10, 22, 28, 15},
		new int[] {14, 15, 28, 27, 16, 0, 0, 0},
		new int[] {9, 6, 19, 22, 12, 11, 24, 25},
		new int[] {12, 11, 24, 25, 13, 0, 0, 0},
		new int[] {19, 24, 27, 22, 25, 28, 0, 0},
		new int[] {16, 28, 27, 13, 25, 24, 0, 0}
	};
	public static HPRef_Struct reftet_1e_3vb = new HPRef_Struct(HP_TET, reftet_1e_3vb_splitedges, reftet_1e_3vb_splitfaces, 0, reftet_1e_3vb_newelstypes, reftet_1e_3vb_newels);






	/*
	// HP_TET_1E_4V
	int reftet_1e_4v_splitedges[][3] =
	{
	  { 1, 2, 5 },
	  { 1, 3, 6 },
	  { 1, 4, 7 },
	  { 2, 1, 8 },
	  { 2, 3, 9 },
	  { 2, 4, 10 },
	  { 3, 1, 11 },
	  { 3, 2, 12 },
	  { 3, 4, 13 },
	  { 4, 1, 14 },
	  { 4, 2, 15 },
	  { 4, 3, 16 },
	  { 0, 0, 0 }
	};
	int reftet_1e_4v_splitfaces[][4] =
	  {
	    { 1, 2, 3, 17 },
	    { 1, 2, 4, 18 },
	    { 1, 3, 4, 19 },
	
	    { 2, 1, 3, 20 },
	    { 2, 1, 4, 21 },
	    { 2, 3, 4, 22 },
	
	    { 3, 1, 2, 23 },
	    { 3, 1, 4, 24 },
	    { 3, 2, 4, 25 },
	
	    { 4, 1, 2, 26 },
	    { 4, 1, 3, 27 },
	    { 4, 2, 3, 28 },
	    { 0, 0, 0, 0 },
	  };
	int reftet_1e_4v_splitelements[][5] =
	  {
	    { 1, 2, 3, 4, 29 },
	    { 2, 3, 4, 1, 30 },
	    { 3, 4, 1, 2, 31 },
	    { 4, 1, 2, 3, 32 },
	    { 0 },
	  };
	HPREF_ELEMENT_TYPE reftet_1e_4v_newelstypes[] =
	{
	  HP_HEX_1E_1V,
	  HP_HEX_1E_1V,
	  HP_HEX_0E_1V,
	  HP_HEX_0E_1V,
	  HP_PRISM_SINGEDGE, HP_PRISM, 
	  HP_PRISM, HP_PRISM, 
	  HP_PRISM, HP_PRISM, 
	  HP_PRISM, HP_PRISM, 
	  HP_PRISM, HP_PRISM, 
	  HP_PRISM, HP_PRISM, 
	  HP_PRISM,
	  HP_PRISM,
	  HP_PRISM,
	  HP_PRISM,
	  HP_TET,
	  HP_NONE,
	};
	int reftet_1e_4v_newels[][8] =
	{
	  { 1, 5, 17, 6, 7, 18, 29, 19 },
	  //  { 2, 9, 20, 8, 10, 22, 30, 21 },
	  { 2, 8, 21, 10, 9, 20, 30, 22 },
	  { 3, 11, 23, 12, 13, 24, 31, 25 },
	  { 4, 15, 26, 14, 16, 28, 32, 27 },
	  { 5, 17, 18, 8, 20, 21 },
	  { 18, 17, 29, 21, 20, 30 },
	  { 6, 19, 17,  11, 24, 23 },
	  { 17, 19, 29,  23, 24, 31 },
	  { 7, 18, 19, 14, 26, 27 },
	  { 19, 18, 29, 27, 26, 32 },
	  { 9, 20, 22, 12, 23, 25 },
	  { 22, 20, 30, 25, 23, 31 },
	  { 10, 22, 21, 15, 28, 26 },
	  { 21, 22, 30, 26, 28, 32 },
	  { 13, 24, 25, 16, 27, 28 },
	  { 25, 24, 31, 28, 27, 32 },
	  { 17, 20, 23, 29, 30, 31 },
	  { 18, 26, 21, 29, 32, 30 },
	  { 19, 24, 27, 29, 31, 32 },
	  { 22, 28, 25, 30, 32, 31 },
	
	  { 29, 30, 31, 32 },
	};
	HPRef_Struct reftet_1e_4v =
	{
	  HP_TET,
	  reftet_1e_4v_splitedges, 
	  reftet_1e_4v_splitfaces, 
	  reftet_1e_4v_splitelements, 
	  reftet_1e_4v_newelstypes, 
	  reftet_1e_4v_newels
	};
	*/




	// HP_TET_1E_4V
	public static int[][] reftet_1e_4v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_1e_4v_splitfaces =
	{
		new int[] {1, 3, 4, 17},
		new int[] {2, 3, 4, 18},
		new int[] {3, 1, 4, 19},
		new int[] {3, 2, 4, 20},
		new int[] {4, 1, 3, 21},
		new int[] {4, 2, 3, 22},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_1e_4v_newelstypes = {HP_TET_1E_1VA, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_HEX, HP_HEX, HP_PRISM, HP_PRISM, HP_PYRAMID, HP_TET_0E_1V, HP_PYRAMID, HP_TET_0E_1V, HP_NONE};

	public static int[][] reftet_1e_4v_newels =
	{
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {7, 6, 17, 10, 9, 18, 0, 0},
		new int[] {7, 10, 18, 17, 14, 15, 22, 21},
		new int[] {9, 6, 17, 18, 12, 11, 19, 20},
		new int[] {17, 19, 21, 18, 20, 22, 0, 0},
		new int[] {16, 22, 21, 13, 20, 19, 0, 0},
		new int[] {14, 15, 22, 21, 16, 0, 0, 0},
		new int[] {4, 14, 16, 15, 0, 0, 0, 0},
		new int[] {12, 11, 19, 20, 13, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {1, 5, 17, 6, 7, 18, 29, 19},
		new int[] {2, 8, 21, 10, 9, 20, 30, 22},
		new int[] {3, 11, 23, 12, 13, 24, 31, 25},
		new int[] {4, 15, 26, 14, 16, 28, 32, 27},
		new int[] {5, 17, 18, 8, 20, 21, 0, 0},
		new int[] {18, 17, 29, 21, 20, 30, 0, 0},
		new int[] {6, 19, 17, 11, 24, 23, 0, 0},
		new int[] {17, 19, 29, 23, 24, 31, 0, 0},
		new int[] {7, 18, 19, 14, 26, 27, 0, 0},
		new int[] {19, 18, 29, 27, 26, 32, 0, 0},
		new int[] {9, 20, 22, 12, 23, 25, 0, 0},
		new int[] {22, 20, 30, 25, 23, 31, 0, 0},
		new int[] {10, 22, 21, 15, 28, 26, 0, 0},
		new int[] {21, 22, 30, 26, 28, 32, 0, 0},
		new int[] {13, 24, 25, 16, 27, 28, 0, 0},
		new int[] {25, 24, 31, 28, 27, 32, 0, 0},
		new int[] {17, 20, 23, 29, 30, 31, 0, 0},
		new int[] {18, 26, 21, 29, 32, 30, 0, 0},
		new int[] {19, 24, 27, 29, 31, 32, 0, 0},
		new int[] {22, 28, 25, 30, 32, 31, 0, 0},
		new int[] {29, 30, 31, 32, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1e_4v = new HPRef_Struct(HP_TET, reftet_1e_4v_splitedges, reftet_1e_4v_splitfaces, 0, reftet_1e_4v_newelstypes, reftet_1e_4v_newels);













	//  HP_TET_2EA_0V,  // 2 edges connected
	public static int[][] reftet_2ea_0v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_0v_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_0v_newelstypes = {HP_PYRAMID_EDGES, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_0v_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {5, 17, 7, 2, 9, 10, 0, 0},
		new int[] {6, 7, 17, 3, 13, 12, 0, 0},
		new int[] {17, 9, 12, 7, 10, 13, 0, 0},
		new int[] {7, 10, 13, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_0v = new HPRef_Struct(HP_TET, reftet_2ea_0v_splitedges, reftet_2ea_0v_splitfaces, 0, reftet_2ea_0v_newelstypes, reftet_2ea_0v_newels);






	//  HP_TET_2EA_1VA,  // 2 edges connected
	public static int[][] reftet_2ea_1va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_1va_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_1va_newelstypes = {HP_PYRAMID_EDGES, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_1va_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 10, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {6, 7, 17, 3, 13, 12, 0, 0},
		new int[] {17, 9, 12, 7, 10, 13, 0, 0},
		new int[] {7, 10, 13, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_1va = new HPRef_Struct(HP_TET, reftet_2ea_1va_splitedges, reftet_2ea_1va_splitfaces, 0, reftet_2ea_1va_newelstypes, reftet_2ea_1va_newels);








	//  HP_TET_2EA_1VB, 
	public static int[][] reftet_2ea_1vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_1vb_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_1vb_newelstypes = {HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_1vb_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {5, 17, 7, 2, 9, 10, 0, 0},
		new int[] {6, 7, 17, 11, 13, 12, 0, 0},
		new int[] {17, 9, 12, 7, 10, 13, 0, 0},
		new int[] {7, 10, 13, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_1vb = new HPRef_Struct(HP_TET, reftet_2ea_1vb_splitedges, reftet_2ea_1vb_splitfaces, 0, reftet_2ea_1vb_newelstypes, reftet_2ea_1vb_newels);






	//  HP_TET_2EA_1VC,  // 2 edges connected
	public static int[][] reftet_2ea_1vc_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_1vc_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 3, 4, 18},
		new int[] {3, 4, 2, 19},
		new int[] {4, 2, 3, 20},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_2ea_1vc_splitelements =
	{
		new int[] {1, 2, 3, 4, 21},
		new int[] {0, 0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_1vc_newelstypes = {HP_PYRAMID_EDGES, HP_TET_0E_1V, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET, HP_TET, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_1vc_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {4, 15, 14, 16, 0, 0, 0, 0},
		new int[] {5, 17, 7, 2, 9, 10, 0, 0},
		new int[] {6, 7, 17, 3, 13, 12, 0, 0},
		new int[] {9, 10, 18, 21, 0, 0, 0, 0},
		new int[] {13, 12, 19, 21, 0, 0, 0, 0},
		new int[] {15, 16, 20, 21, 0, 0, 0, 0},
		new int[] {18, 20, 19, 21, 0, 0, 0, 0},
		new int[] {10, 15, 20, 18, 21, 0, 0, 0},
		new int[] {13, 19, 20, 16, 21, 0, 0, 0},
		new int[] {9, 18, 19, 12, 21, 0, 0, 0},
		new int[] {7, 13, 16, 14, 21, 0, 0, 0},
		new int[] {7, 14, 15, 10, 21, 0, 0, 0},
		new int[] {9, 12, 17, 21, 0, 0, 0, 0},
		new int[] {7, 10, 9, 17, 21, 0, 0, 0},
		new int[] {7, 17, 12, 13, 21, 0, 0, 0},
		new int[] {14, 16, 15, 21, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_1vc = new HPRef_Struct(HP_TET, reftet_2ea_1vc_splitedges, reftet_2ea_1vc_splitfaces, reftet_2ea_1vc_splitelements, reftet_2ea_1vc_newelstypes, reftet_2ea_1vc_newels);












	//  HP_TET_2EA_2VA, 
	public static int[][] reftet_2ea_2va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_2va_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_2va_newelstypes = {HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_2va_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 10, 0, 0},
		new int[] {6, 7, 17, 11, 13, 12, 0, 0},
		new int[] {17, 9, 12, 7, 10, 13, 0, 0},
		new int[] {7, 10, 13, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_2va = new HPRef_Struct(HP_TET, reftet_2ea_2va_splitedges, reftet_2ea_2va_splitfaces, 0, reftet_2ea_2va_newelstypes, reftet_2ea_2va_newels);











	//  HP_TET_2EA_2VB,  // 2 edges connected
	public static int[][] reftet_2ea_2vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_2vb_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 3, 4, 18},
		new int[] {3, 4, 2, 19},
		new int[] {4, 2, 3, 20},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_2ea_2vb_splitelements =
	{
		new int[] {1, 2, 3, 4, 21},
		new int[] {0, 0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_2vb_newelstypes = {HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_TET_0E_1V, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET, HP_TET, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_2vb_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {4, 15, 14, 16, 0, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 10, 0, 0},
		new int[] {6, 7, 17, 3, 13, 12, 0, 0},
		new int[] {9, 10, 18, 21, 0, 0, 0, 0},
		new int[] {13, 12, 19, 21, 0, 0, 0, 0},
		new int[] {15, 16, 20, 21, 0, 0, 0, 0},
		new int[] {18, 20, 19, 21, 0, 0, 0, 0},
		new int[] {10, 15, 20, 18, 21, 0, 0, 0},
		new int[] {13, 19, 20, 16, 21, 0, 0, 0},
		new int[] {9, 18, 19, 12, 21, 0, 0, 0},
		new int[] {7, 13, 16, 14, 21, 0, 0, 0},
		new int[] {7, 14, 15, 10, 21, 0, 0, 0},
		new int[] {9, 12, 17, 21, 0, 0, 0, 0},
		new int[] {7, 10, 9, 17, 21, 0, 0, 0},
		new int[] {7, 17, 12, 13, 21, 0, 0, 0},
		new int[] {14, 16, 15, 21, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_2vb = new HPRef_Struct(HP_TET, reftet_2ea_2vb_splitedges, reftet_2ea_2vb_splitfaces, reftet_2ea_2vb_splitelements, reftet_2ea_2vb_newelstypes, reftet_2ea_2vb_newels);










	//  HP_TET_2EA_2VC,  // 2 edges connected
	public static int[][] reftet_2ea_2vc_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_2vc_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 3, 4, 18},
		new int[] {3, 4, 2, 19},
		new int[] {4, 2, 3, 20},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_2ea_2vc_splitelements =
	{
		new int[] {1, 2, 3, 4, 21},
		new int[] {0, 0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_2vc_newelstypes = {HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_TET_0E_1V, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET, HP_TET, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_2vc_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {4, 15, 14, 16, 0, 0, 0, 0},
		new int[] {5, 17, 7, 2, 9, 10, 0, 0},
		new int[] {6, 7, 17, 11, 13, 12, 0, 0},
		new int[] {9, 10, 18, 21, 0, 0, 0, 0},
		new int[] {13, 12, 19, 21, 0, 0, 0, 0},
		new int[] {15, 16, 20, 21, 0, 0, 0, 0},
		new int[] {18, 20, 19, 21, 0, 0, 0, 0},
		new int[] {10, 15, 20, 18, 21, 0, 0, 0},
		new int[] {13, 19, 20, 16, 21, 0, 0, 0},
		new int[] {9, 18, 19, 12, 21, 0, 0, 0},
		new int[] {7, 13, 16, 14, 21, 0, 0, 0},
		new int[] {7, 14, 15, 10, 21, 0, 0, 0},
		new int[] {9, 12, 17, 21, 0, 0, 0, 0},
		new int[] {7, 10, 9, 17, 21, 0, 0, 0},
		new int[] {7, 17, 12, 13, 21, 0, 0, 0},
		new int[] {14, 16, 15, 21, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_2vc = new HPRef_Struct(HP_TET, reftet_2ea_2vc_splitedges, reftet_2ea_2vc_splitfaces, reftet_2ea_2vc_splitelements, reftet_2ea_2vc_newelstypes, reftet_2ea_2vc_newels);








	//  HP_TET_2EA_3V,  // 2 edges connected
	public static int[][] reftet_2ea_3v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_2ea_3v_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 3, 4, 18},
		new int[] {3, 4, 2, 19},
		new int[] {4, 2, 3, 20},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_2ea_3v_splitelements =
	{
		new int[] {1, 2, 3, 4, 21},
		new int[] {0, 0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_2ea_3v_newelstypes = {HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_0E_1V, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET, HP_TET, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_NONE};
	public static int[][] reftet_2ea_3v_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {4, 15, 14, 16, 0, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 10, 0, 0},
		new int[] {6, 7, 17, 11, 13, 12, 0, 0},
		new int[] {9, 10, 18, 21, 0, 0, 0, 0},
		new int[] {13, 12, 19, 21, 0, 0, 0, 0},
		new int[] {15, 16, 20, 21, 0, 0, 0, 0},
		new int[] {18, 20, 19, 21, 0, 0, 0, 0},
		new int[] {10, 15, 20, 18, 21, 0, 0, 0},
		new int[] {13, 19, 20, 16, 21, 0, 0, 0},
		new int[] {9, 18, 19, 12, 21, 0, 0, 0},
		new int[] {7, 13, 16, 14, 21, 0, 0, 0},
		new int[] {7, 14, 15, 10, 21, 0, 0, 0},
		new int[] {9, 12, 17, 21, 0, 0, 0, 0},
		new int[] {7, 10, 9, 17, 21, 0, 0, 0},
		new int[] {7, 17, 12, 13, 21, 0, 0, 0},
		new int[] {14, 16, 15, 21, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2ea_3v = new HPRef_Struct(HP_TET, reftet_2ea_3v_splitedges, reftet_2ea_3v_splitfaces, reftet_2ea_3v_splitelements, reftet_2ea_3v_newelstypes, reftet_2ea_3v_newels);







	//  HP_TET_2EB_0V,  // 2 opposite edges
	public static int[][] reftet_2eb_0v_splitedges =
	{
		new int[] {1, 3, 5},
		new int[] {1, 4, 6},
		new int[] {2, 3, 7},
		new int[] {2, 4, 8},
		new int[] {3, 1, 9},
		new int[] {3, 2, 10},
		new int[] {4, 1, 11},
		new int[] {4, 2, 12},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_0v_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_0v_newels =
	{
		new int[] {1, 5, 6, 2, 7, 8, 0, 0},
		new int[] {3, 9, 10, 4, 11, 12, 0, 0},
		new int[] {6, 11, 12, 8, 5, 9, 10, 7}
	};
	public static HPRef_Struct reftet_2eb_0v = new HPRef_Struct(HP_TET, reftet_2eb_0v_splitedges, 0, 0, reftet_2eb_0v_newelstypes, reftet_2eb_0v_newels);


	//  HP_TET_2EB_1V,    // V1


	//  HP_TET_2EB_1V,  // 2 opposite edges, V1
	public static int[][] reftet_2eb_1v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_1v_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_1v_newels =
	{
		new int[] {5, 6, 7, 2, 9, 10, 0, 0},
		new int[] {4, 15, 14, 3, 12, 11, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {7, 14, 15, 10, 6, 11, 12, 9}
	};
	public static HPRef_Struct reftet_2eb_1v = new HPRef_Struct(HP_TET, reftet_2eb_1v_splitedges, 0, 0, reftet_2eb_1v_newelstypes, reftet_2eb_1v_newels);



	//  HP_TET_2EB_2VA,  // 2 opposite edges, V1,2
	public static int[][] reftet_2eb_2va_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_2va_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_2va_newels =
	{
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {4, 15, 14, 3, 12, 11, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {7, 14, 15, 10, 6, 11, 12, 9}
	};
	public static HPRef_Struct reftet_2eb_2va = new HPRef_Struct(HP_TET, reftet_2eb_2va_splitedges, 0, 0, reftet_2eb_2va_newelstypes, reftet_2eb_2va_newels);


	//  HP_TET_2EB_2VB,   // V1,3
	public static int[][] reftet_2eb_2vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_2vb_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_2vb_newels =
	{
		new int[] {5, 6, 7, 2, 9, 10, 0, 0},
		new int[] {4, 15, 14, 13, 12, 11, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {3, 13, 11, 12, 0, 0, 0, 0},
		new int[] {7, 14, 15, 10, 6, 11, 12, 9}
	};
	public static HPRef_Struct reftet_2eb_2vb = new HPRef_Struct(HP_TET, reftet_2eb_2vb_splitedges, 0, 0, reftet_2eb_2vb_newelstypes, reftet_2eb_2vb_newels);




	//  HP_TET_2EB_2VC,   // V1,4
	public static int[][] reftet_2eb_2vc_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_2vc_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_2vc_newels =
	{
		new int[] {5, 6, 7, 2, 9, 10, 0, 0},
		new int[] {16, 15, 14, 3, 12, 11, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {4, 16, 15, 14, 0, 0, 0, 0},
		new int[] {7, 14, 15, 10, 6, 11, 12, 9}
	};
	public static HPRef_Struct reftet_2eb_2vc = new HPRef_Struct(HP_TET, reftet_2eb_2vc_splitedges, 0, 0, reftet_2eb_2vc_newelstypes, reftet_2eb_2vc_newels);






	//  HP_TET_2EB_3V,    // V1,2,3
	public static int[][] reftet_2eb_3v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_3v_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_3v_newels =
	{
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {4, 15, 14, 13, 12, 11, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {3, 13, 11, 12, 0, 0, 0, 0},
		new int[] {7, 14, 15, 10, 6, 11, 12, 9}
	};
	public static HPRef_Struct reftet_2eb_3v = new HPRef_Struct(HP_TET, reftet_2eb_3v_splitedges, 0, 0, reftet_2eb_3v_newelstypes, reftet_2eb_3v_newels);






	//  HP_TET_2EB_4V,  // 2 opposite edges
	public static int[][] reftet_2eb_4v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_2eb_4v_newelstypes = {HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_HEX, HP_NONE};
	public static int[][] reftet_2eb_4v_newels =
	{
		new int[] {5, 6, 7, 8, 9, 10, 0, 0},
		new int[] {16, 15, 14, 13, 12, 11, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0},
		new int[] {2, 8, 10, 9, 0, 0, 0, 0},
		new int[] {3, 13, 11, 12, 0, 0, 0, 0},
		new int[] {4, 16, 15, 14, 0, 0, 0, 0},
		new int[] {7, 14, 15, 10, 6, 11, 12, 9}
	};
	public static HPRef_Struct reftet_2eb_4v = new HPRef_Struct(HP_TET, reftet_2eb_4v_splitedges, 0, 0, reftet_2eb_4v_newelstypes, reftet_2eb_4v_newels);

















	//  HP_TET_3EA_0V,  
	public static int[][] reftet_3ea_0v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {4, 2, 12},
		new int[] {4, 3, 13},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ea_0v_splitfaces =
	{
		new int[] {1, 2, 3, 14},
		new int[] {1, 2, 4, 15},
		new int[] {1, 3, 4, 16},
		new int[] {2, 3, 4, 17},
		new int[] {3, 4, 2, 18},
		new int[] {4, 2, 3, 19},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ea_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ea_0v_newelstypes = {HP_HEX_3E_0V, HP_HEX_1E_0V, HP_HEX_1E_0V, HP_HEX_1E_0V, HP_PRISM, HP_PRISM, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_3ea_0v_newels =
	{
		new int[] {1, 5, 14, 6, 7, 15, 20, 16},
		new int[] {5, 2, 8, 14, 15, 9, 17, 20},
		new int[] {3, 6, 14, 10, 11, 16, 20, 18},
		new int[] {7, 4, 12, 15, 16, 13, 19, 20},
		new int[] {11, 13, 16, 18, 19, 20, 0, 0},
		new int[] {15, 12, 9, 20, 19, 17, 0, 0},
		new int[] {8, 10, 14, 17, 18, 20, 0, 0},
		new int[] {20, 17, 18, 19, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ea_0v = new HPRef_Struct(HP_TET, reftet_3ea_0v_splitedges, reftet_3ea_0v_splitfaces, reftet_3ea_0v_splitelements, reftet_3ea_0v_newelstypes, reftet_3ea_0v_newels);










	//  HP_TET_3EA_1V,  
	public static int[][] reftet_3ea_1v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {4, 2, 12},
		new int[] {4, 3, 13},
		new int[] {2, 1, 21},
		new int[] {3, 1, 22},
		new int[] {4, 1, 23},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ea_1v_splitfaces =
	{
		new int[] {1, 2, 3, 14},
		new int[] {1, 2, 4, 15},
		new int[] {1, 3, 4, 16},
		new int[] {2, 3, 4, 17},
		new int[] {3, 4, 2, 18},
		new int[] {4, 2, 3, 19},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ea_1v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ea_1v_newelstypes = {HP_HEX_3E_0V, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_PRISM_SINGEDGE, HP_PRISM, HP_PRISM_SINGEDGE, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_3ea_1v_newels =
	{
		new int[] {1, 5, 14, 6, 7, 15, 20, 16},
		new int[] {2, 21, 9, 8, 0, 0, 0, 0},
		new int[] {5, 14, 15, 21, 8, 9, 0, 0},
		new int[] {15, 14, 20, 9, 8, 17, 0, 0},
		new int[] {6, 16, 14, 3, 11, 10, 0, 0},
		new int[] {14, 16, 20, 10, 11, 18, 0, 0},
		new int[] {7, 15, 16, 4, 12, 13, 0, 0},
		new int[] {16, 15, 20, 13, 12, 19, 0, 0},
		new int[] {11, 13, 16, 18, 19, 20, 0, 0},
		new int[] {15, 12, 9, 20, 19, 17, 0, 0},
		new int[] {8, 10, 14, 17, 18, 20, 0, 0},
		new int[] {20, 17, 18, 19, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ea_1v = new HPRef_Struct(HP_TET, reftet_3ea_1v_splitedges, reftet_3ea_1v_splitfaces, reftet_3ea_1v_splitelements, reftet_3ea_1v_newelstypes, reftet_3ea_1v_newels);










	//  HP_TET_3EA_2V,  
	public static int[][] reftet_3ea_2v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {4, 2, 12},
		new int[] {4, 3, 13},
		new int[] {2, 1, 21},
		new int[] {3, 1, 22},
		new int[] {4, 1, 23},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ea_2v_splitfaces =
	{
		new int[] {1, 2, 3, 14},
		new int[] {1, 2, 4, 15},
		new int[] {1, 3, 4, 16},
		new int[] {2, 3, 4, 17},
		new int[] {3, 4, 2, 18},
		new int[] {4, 2, 3, 19},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ea_2v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ea_2v_newelstypes = {HP_HEX_3E_0V, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_PRISM_SINGEDGE, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_3ea_2v_newels =
	{
		new int[] {1, 5, 14, 6, 7, 15, 20, 16},
		new int[] {2, 21, 9, 8, 0, 0, 0, 0},
		new int[] {5, 14, 15, 21, 8, 9, 0, 0},
		new int[] {15, 14, 20, 9, 8, 17, 0, 0},
		new int[] {3, 22, 10, 11, 0, 0, 0, 0},
		new int[] {6, 16, 14, 22, 11, 10, 0, 0},
		new int[] {14, 16, 20, 10, 11, 18, 0, 0},
		new int[] {7, 15, 16, 4, 12, 13, 0, 0},
		new int[] {16, 15, 20, 13, 12, 19, 0, 0},
		new int[] {11, 13, 16, 18, 19, 20, 0, 0},
		new int[] {15, 12, 9, 20, 19, 17, 0, 0},
		new int[] {8, 10, 14, 17, 18, 20, 0, 0},
		new int[] {20, 17, 18, 19, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ea_2v = new HPRef_Struct(HP_TET, reftet_3ea_2v_splitedges, reftet_3ea_2v_splitfaces, reftet_3ea_2v_splitelements, reftet_3ea_2v_newelstypes, reftet_3ea_2v_newels);








	//  HP_TET_3EA_3V,  
	public static int[][] reftet_3ea_3v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {3, 2, 10},
		new int[] {3, 4, 11},
		new int[] {4, 2, 12},
		new int[] {4, 3, 13},
		new int[] {2, 1, 21},
		new int[] {3, 1, 22},
		new int[] {4, 1, 23},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ea_3v_splitfaces =
	{
		new int[] {1, 2, 3, 14},
		new int[] {1, 2, 4, 15},
		new int[] {1, 3, 4, 16},
		new int[] {2, 3, 4, 17},
		new int[] {3, 4, 2, 18},
		new int[] {4, 2, 3, 19},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ea_3v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ea_3v_newelstypes = {HP_HEX_3E_0V, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_PRISM, HP_PRISM, HP_PRISM, HP_TET, HP_NONE};
	public static int[][] reftet_3ea_3v_newels =
	{
		new int[] {1, 5, 14, 6, 7, 15, 20, 16},
		new int[] {2, 21, 9, 8, 0, 0, 0, 0},
		new int[] {5, 14, 15, 21, 8, 9, 0, 0},
		new int[] {15, 14, 20, 9, 8, 17, 0, 0},
		new int[] {3, 22, 10, 11, 0, 0, 0, 0},
		new int[] {6, 16, 14, 22, 11, 10, 0, 0},
		new int[] {14, 16, 20, 10, 11, 18, 0, 0},
		new int[] {4, 23, 13, 12, 0, 0, 0, 0},
		new int[] {7, 15, 16, 23, 12, 13, 0, 0},
		new int[] {16, 15, 20, 13, 12, 19, 0, 0},
		new int[] {11, 13, 16, 18, 19, 20, 0, 0},
		new int[] {15, 12, 9, 20, 19, 17, 0, 0},
		new int[] {8, 10, 14, 17, 18, 20, 0, 0},
		new int[] {20, 17, 18, 19, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ea_3v = new HPRef_Struct(HP_TET, reftet_3ea_3v_splitedges, reftet_3ea_3v_splitfaces, reftet_3ea_3v_splitelements, reftet_3ea_3v_newelstypes, reftet_3ea_3v_newels);







	//  HP_TET_3EV_0V,  
	public static int[][] reftet_3eb_0v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 4, 13},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3eb_0v_splitfaces =
	{
		new int[] {1, 2, 4, 17},
		new int[] {2, 1, 3, 18},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3eb_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3eb_0v_newelstypes = {HP_PYRAMID_EDGES, HP_PYRAMID_EDGES, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_NONE};
	public static int[][] reftet_3eb_0v_newels =
	{
		new int[] {1, 7, 17, 5, 6, 0, 0, 0},
		new int[] {2, 9, 18, 8, 10, 0, 0, 0},
		new int[] {5, 6, 17, 8, 18, 10, 0, 0},
		new int[] {7, 17, 6, 4, 15, 16, 0, 0},
		new int[] {9, 18, 10, 3, 11, 13, 0, 0},
		new int[] {10, 15, 16, 13, 20, 0, 0, 0},
		new int[] {6, 11, 13, 16, 20, 0, 0, 0},
		new int[] {10, 17, 15, 20, 0, 0, 0, 0},
		new int[] {6, 18, 11, 20, 0, 0, 0, 0},
		new int[] {6, 17, 10, 18, 20, 0, 0, 0},
		new int[] {6, 16, 15, 17, 20, 0, 0, 0},
		new int[] {18, 10, 13, 11, 20, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3eb_0v = new HPRef_Struct(HP_TET, reftet_3eb_0v_splitedges, reftet_3eb_0v_splitfaces, reftet_3eb_0v_splitelements, reftet_3eb_0v_newelstypes, reftet_3eb_0v_newels);









	//  HP_TET_3EV_1V,  
	public static int[][] reftet_3eb_1v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3eb_1v_splitfaces =
	{
		new int[] {1, 2, 4, 17},
		new int[] {2, 1, 3, 18},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3eb_1v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3eb_1v_newelstypes = {HP_PYRAMID_EDGES, HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_NONE};
	public static int[][] reftet_3eb_1v_newels =
	{
		new int[] {1, 7, 17, 5, 6, 0, 0, 0},
		new int[] {2, 9, 18, 8, 10, 0, 0, 0},
		new int[] {3, 12, 13, 11, 0, 0, 0, 0},
		new int[] {5, 6, 17, 8, 18, 10, 0, 0},
		new int[] {7, 17, 6, 4, 15, 16, 0, 0},
		new int[] {9, 18, 10, 12, 11, 13, 0, 0},
		new int[] {10, 15, 16, 13, 20, 0, 0, 0},
		new int[] {6, 11, 13, 16, 20, 0, 0, 0},
		new int[] {10, 17, 15, 20, 0, 0, 0, 0},
		new int[] {6, 18, 11, 20, 0, 0, 0, 0},
		new int[] {6, 17, 10, 18, 20, 0, 0, 0},
		new int[] {6, 16, 15, 17, 20, 0, 0, 0},
		new int[] {18, 10, 13, 11, 20, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3eb_1v = new HPRef_Struct(HP_TET, reftet_3eb_1v_splitedges, reftet_3eb_1v_splitfaces, reftet_3eb_1v_splitelements, reftet_3eb_1v_newelstypes, reftet_3eb_1v_newels);








	//  HP_TET_3EV_2V,  
	public static int[][] reftet_3eb_2v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3eb_2v_splitfaces =
	{
		new int[] {1, 2, 4, 17},
		new int[] {2, 1, 3, 18},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3eb_2v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3eb_2v_newelstypes = {HP_PYRAMID_EDGES, HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_NONE};
	public static int[][] reftet_3eb_2v_newels =
	{
		new int[] {1, 7, 17, 5, 6, 0, 0, 0},
		new int[] {2, 9, 18, 8, 10, 0, 0, 0},
		new int[] {3, 12, 13, 11, 0, 0, 0, 0},
		new int[] {4, 14, 16, 15, 0, 0, 0, 0},
		new int[] {5, 6, 17, 8, 18, 10, 0, 0},
		new int[] {7, 17, 6, 14, 15, 16, 0, 0},
		new int[] {9, 18, 10, 12, 11, 13, 0, 0},
		new int[] {10, 15, 16, 13, 20, 0, 0, 0},
		new int[] {6, 11, 13, 16, 20, 0, 0, 0},
		new int[] {10, 17, 15, 20, 0, 0, 0, 0},
		new int[] {6, 18, 11, 20, 0, 0, 0, 0},
		new int[] {6, 17, 10, 18, 20, 0, 0, 0},
		new int[] {6, 16, 15, 17, 20, 0, 0, 0},
		new int[] {18, 10, 13, 11, 20, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3eb_2v = new HPRef_Struct(HP_TET, reftet_3eb_2v_splitedges, reftet_3eb_2v_splitfaces, reftet_3eb_2v_splitelements, reftet_3eb_2v_newelstypes, reftet_3eb_2v_newels);













	//  HP_TET_3EC_0V,  
	public static int[][] reftet_3ec_0v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ec_0v_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 1, 4, 18},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ec_0v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ec_0v_newelstypes = {HP_PYRAMID_EDGES, HP_PYRAMID_EDGES, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_NONE};
	public static int[][] reftet_3ec_0v_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {2, 8, 18, 10, 9, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 18, 0, 0},
		new int[] {6, 7, 17, 3, 13, 12, 0, 0},
		new int[] {10, 9, 18, 4, 16, 14, 0, 0},
		new int[] {9, 16, 13, 12, 20, 0, 0, 0},
		new int[] {7, 13, 16, 14, 20, 0, 0, 0},
		new int[] {7, 14, 18, 20, 0, 0, 0, 0},
		new int[] {9, 12, 17, 20, 0, 0, 0, 0},
		new int[] {17, 7, 18, 9, 20, 0, 0, 0},
		new int[] {7, 17, 12, 13, 20, 0, 0, 0},
		new int[] {9, 18, 14, 16, 20, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ec_0v = new HPRef_Struct(HP_TET, reftet_3ec_0v_splitedges, reftet_3ec_0v_splitfaces, reftet_3ec_0v_splitelements, reftet_3ec_0v_newelstypes, reftet_3ec_0v_newels);









	//  HP_TET_3EC_1V,  
	public static int[][] reftet_3ec_1v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ec_1v_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 1, 4, 18},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ec_1v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ec_1v_newelstypes = {HP_PYRAMID_EDGES, HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_NONE};
	public static int[][] reftet_3ec_1v_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {2, 8, 18, 10, 9, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 18, 0, 0},
		new int[] {6, 7, 17, 11, 13, 12, 0, 0},
		new int[] {10, 9, 18, 4, 16, 14, 0, 0},
		new int[] {9, 16, 13, 12, 20, 0, 0, 0},
		new int[] {7, 13, 16, 14, 20, 0, 0, 0},
		new int[] {7, 14, 18, 20, 0, 0, 0, 0},
		new int[] {9, 12, 17, 20, 0, 0, 0, 0},
		new int[] {17, 7, 18, 9, 20, 0, 0, 0},
		new int[] {7, 17, 12, 13, 20, 0, 0, 0},
		new int[] {9, 18, 14, 16, 20, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ec_1v = new HPRef_Struct(HP_TET, reftet_3ec_1v_splitedges, reftet_3ec_1v_splitfaces, reftet_3ec_1v_splitelements, reftet_3ec_1v_newelstypes, reftet_3ec_1v_newels);








	//  HP_TET_3EC_2V,  
	public static int[][] reftet_3ec_2v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {2, 3, 9},
		new int[] {2, 4, 10},
		new int[] {3, 1, 11},
		new int[] {3, 2, 12},
		new int[] {3, 4, 13},
		new int[] {4, 1, 14},
		new int[] {4, 2, 15},
		new int[] {4, 3, 16},
		new int[] {0, 0, 0}
	};
	public static int[][] reftet_3ec_2v_splitfaces =
	{
		new int[] {1, 2, 3, 17},
		new int[] {2, 1, 4, 18},
		new int[] {0, 0, 0, 0}
	};
	public static int[][] reftet_3ec_2v_splitelements =
	{
		new int[] {1, 2, 3, 4, 20},
		new int[] {0, 0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftet_3ec_2v_newelstypes = {HP_PYRAMID_EDGES, HP_PYRAMID_EDGES, HP_TET_1E_1VA, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PRISM_SINGEDGE, HP_PYRAMID, HP_PYRAMID, HP_TET, HP_TET, HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, HP_NONE};
	public static int[][] reftet_3ec_2v_newels =
	{
		new int[] {1, 5, 17, 6, 7, 0, 0, 0},
		new int[] {2, 8, 18, 10, 9, 0, 0, 0},
		new int[] {3, 11, 12, 13, 0, 0, 0, 0},
		new int[] {4, 15, 14, 16, 0, 0, 0, 0},
		new int[] {5, 17, 7, 8, 9, 18, 0, 0},
		new int[] {6, 7, 17, 11, 13, 12, 0, 0},
		new int[] {10, 9, 18, 15, 16, 14, 0, 0},
		new int[] {9, 16, 13, 12, 20, 0, 0, 0},
		new int[] {7, 13, 16, 14, 20, 0, 0, 0},
		new int[] {7, 14, 18, 20, 0, 0, 0, 0},
		new int[] {9, 12, 17, 20, 0, 0, 0, 0},
		new int[] {17, 7, 18, 9, 20, 0, 0, 0},
		new int[] {7, 17, 12, 13, 20, 0, 0, 0},
		new int[] {9, 18, 14, 16, 20, 0, 0, 0}
	};
	public static HPRef_Struct reftet_3ec_2v = new HPRef_Struct(HP_TET, reftet_3ec_2v_splitedges, reftet_3ec_2v_splitfaces, reftet_3ec_2v_splitelements, reftet_3ec_2v_newelstypes, reftet_3ec_2v_newels);










	/* ************************ 1 singular face ******************** */


	// HP_TET_1F_0E_0V
	public static int[][] reftet_1f_0e_0v_splitedges =
	{
		new int[] {2, 1, 5},
		new int[] {3, 1, 6},
		new int[] {4, 1, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1f_0e_0v_newelstypes = {HP_PRISM_1FA_0E_0V, HP_TET, HP_NONE};
	public static int[][] reftet_1f_0e_0v_newels =
	{
		new int[] {3, 2, 4, 6, 5, 7, 0, 0},
		new int[] {5, 7, 6, 1, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1f_0e_0v = new HPRef_Struct(HP_TET, reftet_1f_0e_0v_splitedges, 0, 0, reftet_1f_0e_0v_newelstypes, reftet_1f_0e_0v_newels);





	// HP_TET_1F_0E_1VA    ... singular vertex in face
	public static int[][] reftet_1f_0e_1va_splitedges =
	{
		new int[] {2, 1, 5},
		new int[] {2, 3, 6},
		new int[] {2, 4, 7},
		new int[] {3, 1, 8},
		new int[] {4, 1, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1f_0e_1va_newelstypes = {HP_HEX_1F_0E_0V, HP_TET_1F_0E_1VA, HP_TET, HP_NONE};
	public static int[][] reftet_1f_0e_1va_newels =
	{
		new int[] {3, 6, 7, 4, 8, 5, 5, 9},
		new int[] {5, 2, 6, 7, 0, 0, 0, 0},
		new int[] {5, 9, 8, 1, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1f_0e_1va = new HPRef_Struct(HP_TET, reftet_1f_0e_1va_splitedges, 0, 0, reftet_1f_0e_1va_newelstypes, reftet_1f_0e_1va_newels);





	// HP_TET_1F_0E_1VB    ... singular vertex not in face
	public static int[][] reftet_1f_0e_1vb_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {1, 3, 6},
		new int[] {1, 4, 7},
		new int[] {2, 1, 8},
		new int[] {3, 1, 9},
		new int[] {4, 1, 10},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftet_1f_0e_1vb_newelstypes = {HP_PRISM_1FA_0E_0V, HP_PRISM, HP_TET_0E_1V, HP_NONE};
	public static int[][] reftet_1f_0e_1vb_newels =
	{
		new int[] {2, 4, 3, 8, 10, 9, 0, 0},
		new int[] {8, 10, 9, 5, 7, 6, 0, 0},
		new int[] {1, 5, 6, 7, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1f_0e_1vb = new HPRef_Struct(HP_TET, reftet_1f_0e_1vb_splitedges, 0, 0, reftet_1f_0e_1vb_newelstypes, reftet_1f_0e_1vb_newels);








	// HP_TET_1F_1EA_0V  ... sing edge is 1..2
	public static int[][] reftet_1f_1ea_0v_splitedges =
	{
		new int[] {1, 3, 5},
		new int[] {1, 4, 6},
		new int[] {2, 1, 7},
		new int[] {2, 3, 8},
		new int[] {2, 4, 9},
		new int[] {3, 1, 10},
		new int[] {4, 1, 11},
		new int[] {0, 0, 0}
	};

	public static int[][] reftet_1f_1ea_0v_splitfaces =
	{
		new int[] {2, 1, 3, 12},
		new int[] {2, 1, 4, 13},
		new int[] {0, 0, 0, 0}
	};


	public static HPREF_ELEMENT_TYPE[] reftet_1f_1ea_0v_newelstypes = {HP_HEX_1F_0E_0V, HP_PYRAMID_1FB_0E_1VA, HP_TET_1E_1VA, HP_PRISM_SINGEDGE, HP_PRISM, HP_NONE};
	public static int[][] reftet_1f_1ea_0v_newels =
	{
		new int[] {3, 8, 9, 4, 10, 12, 13, 11},
		new int[] {8, 9, 13, 12, 2, 0, 0, 0},
		new int[] {2, 7, 13, 12, 0, 0, 0, 0},
		new int[] {7, 13, 12, 1, 6, 5, 0, 0},
		new int[] {6, 11, 13, 5, 10, 12, 0, 0}
	};
	public static HPRef_Struct reftet_1f_1ea_0v = new HPRef_Struct(HP_TET, reftet_1f_1ea_0v_splitedges, reftet_1f_1ea_0v_splitfaces, 0, reftet_1f_1ea_0v_newelstypes, reftet_1f_1ea_0v_newels);








	// HP_TET_1F_1EB_0V     singular edge in face, edge is 2-3
	public static int[][] reftet_1f_1eb_0v_splitedges =
	{
		new int[] {2, 1, 5},
		new int[] {2, 4, 6},
		new int[] {3, 1, 7},
		new int[] {3, 4, 8},
		new int[] {4, 1, 9},
		new int[] {0, 0, 0}
	};


	public static HPREF_ELEMENT_TYPE[] reftet_1f_1eb_0v_newelstypes = {HP_PRISM_1FB_1EA_0V, HP_PRISM_1FA_0E_0V, HP_TET, HP_NONE};
	public static int[][] reftet_1f_1eb_0v_newels =
	{
		new int[] {3, 8, 7, 2, 6, 5, 0, 0},
		new int[] {6, 4, 8, 5, 9, 7, 0, 0},
		new int[] {5, 9, 7, 1, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_1f_1eb_0v = new HPRef_Struct(HP_TET, reftet_1f_1eb_0v_splitedges, 0, 0, reftet_1f_1eb_0v_newelstypes, reftet_1f_1eb_0v_newels);










	/* ************************ 2 singular faces ******************** */


	// HP_TET_2F_0E_0V
	public static int[][] reftet_2f_0e_0v_splitedges =
	{
		new int[] {1, 2, 5},
		new int[] {2, 1, 6},
		new int[] {3, 1, 7},
		new int[] {3, 2, 8},
		new int[] {4, 1, 9},
		new int[] {4, 2, 10},
		new int[] {0, 0, 0}
	};

	public static int[][] reftet_2f_0e_0v_splitfaces =
	{
		new int[] {3, 1, 2, 11},
		new int[] {4, 1, 2, 12},
		new int[] {0, 0, 0, 0}
	};


	public static HPREF_ELEMENT_TYPE[] reftet_2f_0e_0v_newelstypes = {HP_PRISM_1FA_0E_0V, HP_PRISM_1FA_0E_0V, HP_PRISM_1FB_1EA_0V, HP_PRISM_1FB_1EA_0V, HP_TET, HP_NONE};
	public static int[][] reftet_2f_0e_0v_newels =
	{
		new int[] {2, 10, 8, 6, 12, 11, 0, 0},
		new int[] {1, 7, 9, 5, 11, 12, 0, 0},
		new int[] {4, 10, 12, 3, 8, 11, 0, 0},
		new int[] {3, 7, 11, 4, 9, 12, 0, 0},
		new int[] {5, 6, 11, 12, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftet_2f_0e_0v = new HPRef_Struct(HP_TET, reftet_2f_0e_0v_splitedges, reftet_2f_0e_0v_splitfaces, 0, reftet_2f_0e_0v_newelstypes, reftet_2f_0e_0v_newels);


	// HP_TRIG
	public static int[][] reftrig_splitedges =
	{
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_newelstypes = {HP_TRIG, HP_NONE};
	public static int[][] reftrig_newels =
	{
		new int[] {1, 2, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig = new HPRef_Struct(HP_TRIG, reftrig_splitedges, 0, 0, reftrig_newelstypes, reftrig_newels);



	// HP_TRIG_SINGCORNER
	public static int[][] reftrig_singcorner_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singcorner_newelstypes = {HP_TRIG_SINGCORNER, HP_QUAD, HP_NONE};
	public static int[][] reftrig_singcorner_newels =
	{
		new int[] {1, 4, 5, 0, 0, 0, 0, 0},
		new int[] {2, 3, 5, 4, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singcorner = new HPRef_Struct(HP_TRIG, reftrig_singcorner_splitedges, 0, 0, reftrig_singcorner_newelstypes, reftrig_singcorner_newels);


	/*
	// HP_TRIG_SINGCORNER, trigs only
	int reftrig_singcorner_splitedges[][3] =
	{
	  { 1, 2, 4 },
	  { 1, 3, 5 },
	  { 0, 0, 0 }
	};
	HPREF_ELEMENT_TYPE reftrig_singcorner_newelstypes[] =
	{
	  HP_TRIG_SINGCORNER,
	  HP_TRIG,
	  HP_TRIG,
	  HP_NONE,
	};
	int reftrig_singcorner_newels[][8] =
	{
	  { 1, 4, 5 },
	  { 4, 2, 5 },
	  { 5, 2, 3 },
	};
	HPRef_Struct reftrig_singcorner =
	{
	  HP_TRIG,
	  reftrig_singcorner_splitedges, 
	  0, 0,
	  reftrig_singcorner_newelstypes, 
	  reftrig_singcorner_newels
	};
	*/





	// HP_TRIG_SINGCORNER12
	public static int[][] reftrig_singcorner12_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singcorner12_newelstypes = {HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singcorner12_newels =
	{
		new int[] {1, 4, 5, 0, 0, 0, 0, 0},
		new int[] {2, 7, 6, 0, 0, 0, 0, 0},
		new int[] {4, 6, 7, 5, 0, 0, 0, 0},
		new int[] {5, 7, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singcorner12 = new HPRef_Struct(HP_TRIG, reftrig_singcorner12_splitedges, 0, 0, reftrig_singcorner12_newelstypes, reftrig_singcorner12_newels);




	// HP_TRIG_SINGCORNER123_2D
	public static int[][] reftrig_singcorner123_2D_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {3, 1, 8},
		new int[] {3, 2, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singcorner123_2D_newelstypes = {HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_TRIG_SINGCORNER, HP_QUAD, HP_QUAD, HP_NONE};
	public static int[][] reftrig_singcorner123_2D_newels =
	{
		new int[] {1, 4, 5, 0, 0, 0, 0, 0},
		new int[] {2, 7, 6, 0, 0, 0, 0, 0},
		new int[] {3, 8, 9, 0, 0, 0, 0, 0},
		new int[] {4, 6, 8, 5, 0, 0, 0, 0},
		new int[] {6, 7, 9, 8, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singcorner123_2D = new HPRef_Struct(HP_TRIG, reftrig_singcorner123_2D_splitedges, 0, 0, reftrig_singcorner123_2D_newelstypes, reftrig_singcorner123_2D_newels);






	// HP_TRIG_SINGCORNER123
	public static int[][] reftrig_singcorner123_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {3, 1, 8},
		new int[] {3, 2, 9},
		new int[] {0, 0, 0}
	};

	public static int[][] reftrig_singcorner123_splitfaces =
	{
		new int[] {1, 2, 3, 10},
		new int[] {2, 3, 1, 11},
		new int[] {3, 1, 2, 12},
		new int[] {0, 0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singcorner123_newelstypes = {HP_DUMMY_QUAD_SINGCORNER, HP_DUMMY_QUAD_SINGCORNER, HP_DUMMY_QUAD_SINGCORNER, HP_QUAD, HP_QUAD, HP_QUAD, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singcorner123_newels =
	{
		new int[] {1, 4, 10, 5, 0, 0, 0, 0},
		new int[] {2, 7, 11, 6, 0, 0, 0, 0},
		new int[] {3, 8, 12, 9, 0, 0, 0, 0},
		new int[] {4, 6, 11, 10, 0, 0, 0, 0},
		new int[] {7, 9, 12, 11, 0, 0, 0, 0},
		new int[] {8, 5, 10, 12, 0, 0, 0, 0},
		new int[] {10, 11, 12, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singcorner123 = new HPRef_Struct(HP_TRIG, reftrig_singcorner123_splitedges, reftrig_singcorner123_splitfaces, 0, reftrig_singcorner123_newelstypes, reftrig_singcorner123_newels);

	// HP_TRIG_SINGEDGE
	public static int[][] reftrig_singedge_splitedges =
	{
		new int[] {2, 3, 4},
		new int[] {1, 3, 5},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedge_newelstypes = {HP_TRIG, HP_QUAD_SINGEDGE, HP_NONE};
	public static int[][] reftrig_singedge_newels =
	{
		new int[] {4, 3, 5, 0, 0, 0, 0, 0},
		new int[] {1, 2, 4, 5, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedge = new HPRef_Struct(HP_TRIG, reftrig_singedge_splitedges, 0, 0, reftrig_singedge_newelstypes, reftrig_singedge_newels);






	// HP_TRIG_SINGEDGECORNER1
	public static int[][] reftrig_singedgecorner1_splitedges =
	{
		new int[] {1, 2, 6},
		new int[] {1, 3, 5},
		new int[] {2, 3, 4},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner1_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_QUAD_SINGEDGE, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedgecorner1_newels =
	{
		new int[] {1, 6, 5, 0, 0, 0, 0, 0},
		new int[] {6, 2, 4, 5, 0, 0, 0, 0},
		new int[] {5, 4, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner1 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner1_splitedges, 0, 0, reftrig_singedgecorner1_newelstypes, reftrig_singedgecorner1_newels);








	// HP_TRIG_SINGEDGECORNER2
	public static int[][] reftrig_singedgecorner2_splitedges =
	{
		new int[] {2, 1, 6},
		new int[] {1, 3, 5},
		new int[] {2, 3, 4},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner2_newelstypes = {HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedgecorner2_newels =
	{
		new int[] {6, 2, 4, 0, 0, 0, 0, 0},
		new int[] {1, 6, 4, 5, 0, 0, 0, 0},
		new int[] {5, 4, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner2 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner2_splitedges, 0, 0, reftrig_singedgecorner2_newelstypes, reftrig_singedgecorner2_newels);




	// HP_TRIG_SINGEDGECORNER12
	public static int[][] reftrig_singedgecorner12_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner12_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedgecorner12_newels =
	{
		new int[] {1, 4, 5, 0, 0, 0, 0, 0},
		new int[] {6, 2, 7, 0, 0, 0, 0, 0},
		new int[] {4, 6, 7, 5, 0, 0, 0, 0},
		new int[] {5, 7, 3, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner12 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner12_splitedges, 0, 0, reftrig_singedgecorner12_newelstypes, reftrig_singedgecorner12_newels);







	// HP_TRIG_SINGEDGECORNER3
	public static int[][] reftrig_singedgecorner3_splitedges =
	{
		new int[] {1, 3, 4},
		new int[] {3, 1, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner3_newelstypes = {HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] reftrig_singedgecorner3_newels =
	{
		new int[] {1, 2, 6, 4, 0, 0, 0, 0},
		new int[] {4, 6, 7, 5, 0, 0, 0, 0},
		new int[] {3, 5, 7, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner3 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner3_splitedges, 0, 0, reftrig_singedgecorner3_newelstypes, reftrig_singedgecorner3_newels);




	// HP_TRIG_SINGEDGECORNER13
	public static int[][] reftrig_singedgecorner13_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 3, 6},
		new int[] {3, 1, 7},
		new int[] {3, 2, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner13_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] reftrig_singedgecorner13_newels =
	{
		new int[] {1, 4, 5, 0, 0, 0, 0, 0},
		new int[] {4, 2, 6, 5, 0, 0, 0, 0},
		new int[] {5, 6, 8, 7, 0, 0, 0, 0},
		new int[] {3, 7, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner13 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner13_splitedges, 0, 0, reftrig_singedgecorner13_newelstypes, reftrig_singedgecorner13_newels);





	// HP_TRIG_SINGEDGECORNER23
	public static int[][] reftrig_singedgecorner23_splitedges =
	{
		new int[] {1, 3, 4},
		new int[] {2, 1, 5},
		new int[] {2, 3, 6},
		new int[] {3, 1, 7},
		new int[] {3, 2, 8},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner23_newelstypes = {HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] reftrig_singedgecorner23_newels =
	{
		new int[] {5, 2, 6, 0, 0, 0, 0, 0},
		new int[] {1, 5, 6, 4, 0, 0, 0, 0},
		new int[] {4, 6, 8, 7, 0, 0, 0, 0},
		new int[] {3, 7, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner23 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner23_splitedges, 0, 0, reftrig_singedgecorner23_newelstypes, reftrig_singedgecorner23_newels);



	// HP_TRIG_SINGEDGECORNER123
	public static int[][] reftrig_singedgecorner123_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {3, 1, 8},
		new int[] {3, 2, 9},
		new int[] {0, 0, 0}
	};
	public static HPREF_ELEMENT_TYPE[] reftrig_singedgecorner123_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD, HP_TRIG_SINGCORNER, HP_NONE};
	public static int[][] reftrig_singedgecorner123_newels =
	{
		new int[] {1, 4, 5, 0, 0, 0, 0, 0},
		new int[] {6, 2, 7, 0, 0, 0, 0, 0},
		new int[] {4, 6, 7, 5, 0, 0, 0, 0},
		new int[] {5, 7, 9, 8, 0, 0, 0, 0},
		new int[] {3, 8, 9, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedgecorner123 = new HPRef_Struct(HP_TRIG, reftrig_singedgecorner123_splitedges, 0, 0, reftrig_singedgecorner123_newelstypes, reftrig_singedgecorner123_newels);

	// HP_TRIG_SINGEDGES
	public static int[][] reftrig_singedges_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {0, 0, 0}
	};
	public static int[][] reftrig_singedges_splitfaces =
	{
		new int[] {1, 2, 3, 8},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftrig_singedges_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedges_newels =
	{
		new int[] {1, 4, 8, 0, 0, 0, 0, 0},
		new int[] {5, 1, 8, 0, 0, 0, 0, 0},
		new int[] {4, 2, 6, 8, 0, 0, 0, 0},
		new int[] {3, 5, 8, 7, 0, 0, 0, 0},
		new int[] {6, 7, 8, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedges = new HPRef_Struct(HP_TRIG, reftrig_singedges_splitedges, reftrig_singedges_splitfaces, 0, reftrig_singedges_newelstypes, reftrig_singedges_newels);








	// HP_TRIG_SINGEDGES2
	public static int[][] reftrig_singedges2_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {3, 2, 8},
		new int[] {0, 0, 0}
	};
	public static int[][] reftrig_singedges2_splitfaces =
	{
		new int[] {1, 2, 3, 9},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftrig_singedges2_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER2, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedges2_newels =
	{
		new int[] {1, 4, 9, 0, 0, 0, 0, 0},
		new int[] {5, 1, 9, 0, 0, 0, 0, 0},
		new int[] {4, 6, 7, 9, 0, 0, 0, 0},
		new int[] {3, 5, 9, 8, 0, 0, 0, 0},
		new int[] {6, 2, 7, 0, 0, 0, 0, 0},
		new int[] {7, 8, 9, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedges2 = new HPRef_Struct(HP_TRIG, reftrig_singedges2_splitedges, reftrig_singedges2_splitfaces, 0, reftrig_singedges2_newelstypes, reftrig_singedges2_newels);




	// HP_TRIG_SINGEDGES3
	public static int[][] reftrig_singedges3_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 3, 6},
		new int[] {3, 1, 7},
		new int[] {3, 2, 8},
		new int[] {0, 0, 0}
	};
	public static int[][] reftrig_singedges3_splitfaces =
	{
		new int[] {1, 2, 3, 9},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftrig_singedges3_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedges3_newels =
	{
		new int[] {1, 4, 9, 0, 0, 0, 0, 0},
		new int[] {5, 1, 9, 0, 0, 0, 0, 0},
		new int[] {4, 2, 6, 9, 0, 0, 0, 0},
		new int[] {7, 5, 9, 8, 0, 0, 0, 0},
		new int[] {3, 7, 8, 0, 0, 0, 0, 0},
		new int[] {6, 8, 9, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedges3 = new HPRef_Struct(HP_TRIG, reftrig_singedges3_splitedges, reftrig_singedges3_splitfaces, 0, reftrig_singedges3_newelstypes, reftrig_singedges3_newels);






	// HP_TRIG_SINGEDGES23
	public static int[][] reftrig_singedges23_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {1, 3, 5},
		new int[] {2, 1, 6},
		new int[] {2, 3, 7},
		new int[] {3, 1, 8},
		new int[] {3, 2, 9},
		new int[] {0, 0, 0}
	};
	public static int[][] reftrig_singedges23_splitfaces =
	{
		new int[] {1, 2, 3, 10},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftrig_singedges23_newelstypes = {HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_TRIG, HP_NONE};
	public static int[][] reftrig_singedges23_newels =
	{
		new int[] {1, 4, 10, 0, 0, 0, 0, 0},
		new int[] {5, 1, 10, 0, 0, 0, 0, 0},
		new int[] {4, 6, 7, 10, 0, 0, 0, 0},
		new int[] {8, 5, 10, 9, 0, 0, 0, 0},
		new int[] {6, 2, 7, 0, 0, 0, 0, 0},
		new int[] {3, 8, 9, 0, 0, 0, 0, 0},
		new int[] {7, 9, 10, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_singedges23 = new HPRef_Struct(HP_TRIG, reftrig_singedges23_splitedges, reftrig_singedges23_splitfaces, 0, reftrig_singedges23_newelstypes, reftrig_singedges23_newels);


	// HP_TRIG_3SINGEDGES
	public static int[][] reftrig_3singedges_splitedges =
	{
		new int[] {1, 2, 4},
		new int[] {2, 1, 5},
		new int[] {2, 3, 6},
		new int[] {3, 2, 7},
		new int[] {3, 1, 8},
		new int[] {1, 3, 9},
		new int[] {0, 0, 0}
	};
	public static int[][] reftrig_3singedges_splitfaces =
	{
		new int[] {1, 2, 3, 10},
		new int[] {2, 3, 1, 11},
		new int[] {3, 1, 2, 12},
		new int[] {0, 0, 0, 0}
	};

	public static HPREF_ELEMENT_TYPE[] reftrig_3singedges_newelstypes = {HP_TRIG, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_QUAD_SINGEDGE, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_TRIG_SINGEDGECORNER1, HP_TRIG_SINGEDGECORNER2, HP_NONE};
	public static int[][] reftrig_3singedges_newels =
	{
		new int[] {10, 11, 12, 0, 0, 0, 0, 0},
		new int[] {4, 5, 11, 10, 0, 0, 0, 0},
		new int[] {6, 7, 12, 11, 0, 0, 0, 0},
		new int[] {8, 9, 10, 12, 0, 0, 0, 0},
		new int[] {1, 4, 10, 0, 0, 0, 0, 0},
		new int[] {9, 1, 10, 0, 0, 0, 0, 0},
		new int[] {2, 6, 11, 0, 0, 0, 0, 0},
		new int[] {5, 2, 11, 0, 0, 0, 0, 0},
		new int[] {3, 8, 12, 0, 0, 0, 0, 0},
		new int[] {7, 3, 12, 0, 0, 0, 0, 0}
	};
	public static HPRef_Struct reftrig_3singedges = new HPRef_Struct(HP_TRIG, reftrig_3singedges_splitedges, reftrig_3singedges_splitfaces, 0, reftrig_3singedges_newelstypes, reftrig_3singedges_newels);

}

namespace netgen
{
	public static class GlobalMembers
	{
	  // extern DLL_HEADER Array<GeometryRegister*> geometryregister; 
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern DLL_HEADER GeometryRegisterArray geometryregister;

	  public static DLL_HEADER GeometryRegisterArray geometryregister = new DLL_HEADER();

	  internal static RegisterClassForArchive<NetgenGeometry> regnggeo = new RegisterClassForArchive<NetgenGeometry>();

   /*! Philippose - 13/07/2009
       Main function implementing automated assignment of 
       Boundary Condition numbers based on face colours

       This functionality is currently implemtented at the mesh 
       level, and hence allows colour based assignment of boundary 
       conditions for any geometry type within netgen which 
       supports face colours
   */
	   /*! \brief Automatically assign boundary conditions for meshes
	
	       This function allows the boundary condition numbers of a 
	       mesh created in Netgen to be automatically assigned based on 
	       the colours of each face.
	
	       Currently, two algorithms are utilised to assign the BC Properties:
	       1. Automatic assignment using a user defined colour profile file 
	          which defines which RGB colours are to be assigned to which 
	          BC Property number
	          - A default profile file exists in the Netgen folder called 
	            "netgen.ocf"
	       
	       2. The second algorithm uses the following automated algorithm:
	          - Extract all the colours present in the mesh
	          - Use colour index 0 (zero) for all faces with no colour defined
	          - Calculate the number of faces of the surface mesh for each colour
	          - Sort the number of surface elements in ascending order, with the 
	            colour indices as a slave
	          - Use the indices of the sorted array as the BC property number
	
	          Example: If there are 3 colours, present in the mesh and the number 
	          of surface elements for each colour are:
	          - Colour 0: 8500
	          - Colour 1: 120
	          - Colour 2: 2200
	          - Colour 3: 575
	
	          The above is sorted in ascending order and assigned as BC Properties:
	          - BC Prop 0: 120  : Colour 1
	          - BC Prop 1: 575  : Colour 3
	          - BC Prop 2: 2200 : Colour 2
	          - BC Prop 3: 8500 : Colour 0 (no colour defined)
	   */
	   //extern void OCCAutoColourBcProps(Mesh & mesh, OCCGeometry & occgeometry, const char *occcolourfile);
	   public static void AutoColourBcProps(Mesh mesh, string bccolourfile)
	   {
		  // Go directly to the alternate algorithm if no colour profile file was specified
		  if (!bccolourfile)
		  {
			 PrintMessage(1, "AutoColourBcProps: Using Automatic Colour based boundary property assignment algorithm");
			 AutoColourAlg_Sorted(mesh);
		  }
		  else
		  {
			 ifstream ocf = new ifstream(bccolourfile);

			 // If there was an error opening the Colour profile file, jump to the alternate
			 // algorithm after printing a message
			 if (ocf == null)
			 {
				PrintMessage(1, "AutoColourBcProps: Error loading Boundary Colour Profile file ", bccolourfile, " ....", "Switching to Automatic Assignment algorithm!");

				AutoColourAlg_Sorted(mesh);
			 }
			 // If the file opens successfully, call the function which assigns boundary conditions
			 // based on the colour profile file
			 else
			 {
				PrintMessage(1, "AutoColourBcProps: Using Boundary Colour Profile file: ");
				PrintMessage(1, "  ", bccolourfile);
				AutoColourAlg_UserProfile(mesh, ocf);

				// Make sure the file is closed before exiting the function
				if (ocf.is_open())
				{
				   ocf.close();
				}
			 }
		  }
	   }

   /*! Philippose - 11/07/2009
       Function to create a list of all the unique colours 
       available in a given mesh
   */

	   public static void GetFaceColours(Mesh mesh, Array<Vec3d> face_colours)
	   {
		  face_colours.SetSize(1);
		  face_colours.Elem(1) = mesh.GetFaceDescriptor(1).SurfColour();

		  for (int i = 1; i <= mesh.GetNFD(); i++)
		  {
			 Vec3d face_colour = mesh.GetFaceDescriptor(i).SurfColour();
			 bool col_found = false;

			 for (int j = 1; j <= face_colours.Size(); j++)
			 {
				if (ColourMatch(face_colours.Elem(j), new netgen.Vec3d(face_colour)))
				{
				   col_found = true;
				   break;
				}
			 }

			 if (!col_found)
			 {
				 face_colours.Append(face_colour);
			 }
		  }

		  if (printmessage_importance >= 3)
		  {
			 Console.Write("\n");
			 Console.Write("-------- Face Colours --------");
			 Console.Write("\n");
			 for (int i = 1; i <= face_colours.Size(); i++)
			 {
				Console.Write(face_colours.Elem(i));
				Console.Write("\n");
			 }
			 Console.Write("------------------------------");
			 Console.Write("\n");
		  }
	   }

   /*! Philippose - 11/07/2009
       Function to check if two RGB colours are equal

       Note#1: Currently uses unweighted Euclidean Distance 
       for colour matching.

       Note#2: The tolerance used for deciding whether two 
       colours match is defined as "eps" and is currently 
       2.5e-5 (for square of distance)
   */

	   public static bool ColourMatch(Vec3d col1, Vec3d col2, double eps = 2.5e-05)
	   {
		  if (eps <= 0.0)
		  {
			  eps = DefineConstants.DEFAULT_EPS;
		  }

		  bool colmatch = false;

		  if (Dist2(col1, col2) < eps)
		  {
			  colmatch = true;
		  }

		  return colmatch;
	   }
	   // Default colour to be used for boundary condition number "0"

	   // Boundary condition number to use if a face does not have a 
	   // colour assigned to it, or if the colour is the above defined 
	   // default colour

	   // Default tolerance for colour matching (using Euclidean distance)






	   /*! Philippose - 11/07/2009
	       Assign boundary condition numbers based on a user defined 
	       colour profile file.
	
	       The default profile file is "netgen.ocf"
	
	       If the mesh contains colours not defined in the profile,
	       netgen automatically starts assigning each new colour a 
	       new boundary condition number starting from the highest 
	       boundary condition number specified in the profile file.
	   */
	   public static void AutoColourAlg_UserProfile(Mesh mesh, ifstream ocf)
	   {
		  string ocf_inp = new string(new char[100]);
		  bool header_found = false;

		  // Number of colour specifications in the 
		  // user profile file
		  int numentries = 0;
		  while ((ocf.good()) && (!header_found))
		  {
			 ocf >> ocf_inp;
			 if (string.Compare(ocf_inp,"boundary_colours") == 0)
			 {
				 header_found = true;
			 }
		  }

		  if (!header_found)
		  {
			 ocf.close();
			 throw new Exception("AutoColourAlg_UserProfile: Invalid or empty Boundary Colour Profile file\n");
			 return;
		  }

		  // Read in the number of entries from file
		  ocf >> numentries;
		  if (numentries > 0)
		  {
			 if (!ocf.good())
			 {
				ocf.close();
				throw new Exception("AutoColourAlg_UserProfile: Invalid or empty Boundary Colour Profile file\n");
				return;
			 }

			 PrintMessage(3, "Number of colour entries: ", numentries);
		  }
		  else
		  {
			 ocf.close();
			 PrintMessage(3, "AutoColourAlg_UserProfile: No Boundary Colour entries found.... no changes made!");
			 return;
		  }

		  // Arrays to hold the specified RGB colour triplets as well 
		  // as the associated boundary condition number
		  Array<Vec3d> bc_colours = new Array<Vec3d>(numentries);
		  Array<int> bc_num = new Array<int>(numentries);
		  Array<bool> bc_used = new Array<bool>(numentries);

		  // Actually read in the data from the file
		  for (int i = 1; i <= numentries; i++)
		  {
			 int bcnum;
			 // double col_red, col_green, col_blue;

			 ocf >> bcnum;
			 // Boundary condition number DEFAULT_BCNUM is reserved for 
			 // faces which have the default colour Green (0.0,1.0,0.0)
			 // To prevent confusion, no boundary numbery below this default 
			 // are permitted
			 if (bcnum < (DefineConstants.DEFAULT_BCNUM + 1))
			 {
				 bcnum = DefineConstants.DEFAULT_BCNUM + 1;
			 }

			 bc_num.Elem(i) = bcnum;
			 bc_used.Elem(i) = false;
			 ocf >> bc_colours.Elem(i).X() >> bc_colours.Elem(i).Y() >> bc_colours.Elem(i).Z();

			 if (!ocf.good())
			 {
				ocf.close();
				throw new Exception("Boundary Colour file error: Number of entries do not match specified list size!!\n");
				return;
			 }

			 // Bound checking of the values
			 // The RGB values should be between 0.0 and 1.0
			 if (bc_colours.Elem(bcnum).X() < 0.0)
			 {
				 bc_colours.Elem(bcnum).X() = 0.0;
			 }
			 if (bc_colours.Elem(bcnum).X() > 1.0)
			 {
				 bc_colours.Elem(bcnum).X() = 1.0;
			 }
			 if (bc_colours.Elem(bcnum).Y() < 0.0)
			 {
				 bc_colours.Elem(bcnum).X() = 0.0;
			 }
			 if (bc_colours.Elem(bcnum).Y() > 1.0)
			 {
				 bc_colours.Elem(bcnum).X() = 1.0;
			 }
			 if (bc_colours.Elem(bcnum).Z() < 0.0)
			 {
				 bc_colours.Elem(bcnum).X() = 0.0;
			 }
			 if (bc_colours.Elem(bcnum).Z() > 1.0)
			 {
				 bc_colours.Elem(bcnum).X() = 1.0;
			 }
		  }

		  PrintMessage(3, "Successfully loaded Boundary Colour Profile file....");
		  ocf.close();

		  // Find the highest boundary condition number in the list
		  // All colours in the geometry which are not specified in the 
		  // list will be given boundary condition numbers higher than this 
		  // number
		  int max_bcnum = DefineConstants.DEFAULT_BCNUM;
		  for (int i = 1; i <= bc_num.Size();i++)
		  {
			 if (bc_num.Elem(i) > max_bcnum)
			 {
				 max_bcnum = bc_num.Elem(i);
			 }
		  }

		  PrintMessage(3, "Highest boundary number in list = ", max_bcnum);

		  Array<Vec3d> all_colours = new Array<Vec3d>();

		  // Extract all the colours to see how many there are
		  GetFaceColours(mesh, all_colours);
		  PrintMessage(3, "\nNumber of colours defined in Mesh: ", all_colours.Size());

		  if (all_colours.Size() == 0)
		  {
			 PrintMessage(3, "No colour data detected in Mesh... no changes made!");
			 return;
		  }

		  int nfd = mesh.GetNFD();

		  for (int face_index = 1; face_index <= nfd; face_index++)
		  {
			 // Temporary container for individual face colours
			 Vec3d face_colour = new Vec3d();

			 // Get the colour of the face being currently processed
			 face_colour.CopyFrom(mesh.GetFaceDescriptor(face_index).SurfColour());
			 if (!ColourMatch(new netgen.Vec3d(face_colour), new Vec3d(DefineConstants.DEFAULT_R, DefineConstants.DEFAULT_G, DefineConstants.DEFAULT_B)))
			 {
				// Boolean variable to check if the boundary condition was applied 
				// or not... not applied would imply that the colour of the face 
				// does not exist in the list of colours in the profile file
				bool bc_assigned = false;

				for (int col_index = 1; col_index <= bc_colours.Size(); col_index++)
				{
				   if ((ColourMatch(new netgen.Vec3d(face_colour), bc_colours.Elem(col_index))) && (!bc_assigned))
				   {
					  mesh.GetFaceDescriptor(face_index).SetBCProperty(bc_num.Elem(col_index));
					  bc_used.Elem(col_index) = true;
					  bc_assigned = true;
					  break;
				   }
				}

				// If the colour was not found in the list, add it to the list, and assign 
				// the next free boundary condition number to it
				if (!bc_assigned)
				{
				   max_bcnum++;
				   bc_num.Append(max_bcnum);
				   bc_colours.Append(face_colour);
				   bc_used.Append(true);

				   mesh.GetFaceDescriptor(face_index).SetBCProperty(max_bcnum);
				}
			 }
			 else
			 {
				// Set the boundary condition number to the default one
				mesh.GetFaceDescriptor(face_index).SetBCProperty(DefineConstants.DEFAULT_BCNUM);
			 }
		  }

		  // User Information of the results of the operation
		  Vec3d ref_colour = new Vec3d(0.0, 1.0, 0.0);
		  PrintMessage(3, "Colour based Boundary Condition Property details:");
		  for (int bc_index = 0; bc_index <= bc_num.Size(); bc_index++)
		  {
			 if (bc_index > 0)
			 {
				 ref_colour = bc_colours.Elem(bc_index);
			 }

			 if (bc_index == 0)
			 {
				PrintMessage(3, "BC Property: ", DefineConstants.DEFAULT_BCNUM);
				PrintMessage(3, "   RGB Face Colour = ", ref_colour, "", "\n");
			 }
			 else if (bc_used.Elem(bc_index))
			 {
				PrintMessage(3, "BC Property: ", bc_num.Elem(bc_index));
				PrintMessage(3, "   RGB Face Colour = ", ref_colour, "", "\n");
			 }
		  }
	   }





	   /*! Philippose - 11/07/2009
	       Assign boundary condition numbers based on the colours 
	       assigned to each face in the mesh using an automated 
	       algorithm.
	
	       The particular algorithm used has been briefly explained 
	       in the header file "occauxfunctions.hpp"
	   */
	   public static void AutoColourAlg_Sorted(Mesh mesh)
	   {
		  Array<Vec3d> all_colours = new Array<Vec3d>();
		  Array<int> faces_sorted = new Array<int>();
		  Array<int> colours_sorted = new Array<int>();

		  // Extract all the colours to see how many there are
		  GetFaceColours(mesh, all_colours);

		  // Delete the default colour from the list since it will be accounted 
		  // for automatically
		  for (int i = 1; i <= all_colours.Size(); i++)
		  {
			 if (ColourMatch(all_colours.Elem(i), new Vec3d(DefineConstants.DEFAULT_R, DefineConstants.DEFAULT_G, DefineConstants.DEFAULT_B)))
			 {
				all_colours.DeleteElement(i);
				break;
			 }
		  }
		  PrintMessage(3, "\nNumber of colours defined in Mesh: ", all_colours.Size());

		  if (all_colours.Size() == 0)
		  {
			 PrintMessage(3, "No colour data detected in Mesh... no changes made!");
			 return;
		  }

		  // One more slot than the number of colours are required, to 
		  // account for individual faces which have no colour data 
		  // assigned to them in the CAD software
		  faces_sorted.SetSize(all_colours.Size() + 1);
		  colours_sorted.SetSize(all_colours.Size() + 1);
		  faces_sorted = 0;

		  // Slave Array to identify the colours the faces were assigned to, 
		  // after the bubble sort routine to sort the automatic boundary 
		  // identifiers according to the number of surface mesh elements 
		  // of a given colour
		  for (int i = 0; i <= all_colours.Size(); i++)
		  {
			  colours_sorted[i] = i;
		  }

		  // Used to hold the number of surface elements without any OCC 
		  // colour definition
		  int no_colour_faces = 0;

		  // Index in the faces array assigned to faces without any 
		  // or the default colour definition
		  int no_colour_index = 0;

		  int nfd = mesh.GetNFD();

		  // Extract the number of surface elements having a given colour
		  // And save this number into an array for later sorting
		  for (int face_index = 1; face_index <= nfd; face_index++)
		  {
			 Array<SurfaceElementIndex> se_face = new Array<SurfaceElementIndex>();

			 mesh.GetSurfaceElementsOfFace(face_index, se_face);

			 Vec3d face_colour = new Vec3d();

			 face_colour.CopyFrom(mesh.GetFaceDescriptor(face_index).SurfColour());
			 if (!ColourMatch(new netgen.Vec3d(face_colour), new Vec3d(DefineConstants.DEFAULT_R, DefineConstants.DEFAULT_G, DefineConstants.DEFAULT_B)))
			 {
				for (int i = 1; i <= all_colours.Size(); i++)
				{
				   if (ColourMatch(new netgen.Vec3d(face_colour), all_colours.Elem(i)))
				   {
					  faces_sorted[i] = faces_sorted[i] + se_face.Size();
				   }
				}
			 }
			 else
			 {
				// Add the number of surface elements without any colour 
				// definition separately
				no_colour_faces = no_colour_faces + se_face.Size();
			 }
		  }

		  // Sort the face colour indices according to the number of surface 
		  // mesh elements which have a specific colour
		  BubbleSort(new Array<int>(faces_sorted), new Array<int>(colours_sorted));

		  // Now update the array position assigned for surface elements 
		  // without any colour definition with the number of elements
		  faces_sorted[no_colour_index] = no_colour_faces;

		  // Now actually assign the BC Property to the respective faces
		  for (int face_index = 1; face_index <= nfd; face_index++)
		  {
			 Vec3d face_colour = new Vec3d();

			 face_colour.CopyFrom(mesh.GetFaceDescriptor(face_index).SurfColour());
			 if (!ColourMatch(new netgen.Vec3d(face_colour), new Vec3d(DefineConstants.DEFAULT_R, DefineConstants.DEFAULT_G, DefineConstants.DEFAULT_B)))
			 {
				for (int i = 0; i < colours_sorted.Size(); i++)
				{
				   Vec3d ref_colour = new Vec3d();
				   if (i != no_colour_index)
				   {
					   ref_colour = all_colours.Elem(colours_sorted[i]);
				   }

				   if (ColourMatch(new netgen.Vec3d(face_colour), new netgen.Vec3d(ref_colour)))
				   {
					  mesh.GetFaceDescriptor(face_index).SetBCProperty(i + DefineConstants.DEFAULT_BCNUM);
				   }
				}
			 }
			 else
			 {
				mesh.GetFaceDescriptor(face_index).SetBCProperty(DefineConstants.DEFAULT_BCNUM);
			 }

			 PrintMessage(4, "Face number: ", face_index, " ; BC Property = ", mesh.GetFaceDescriptor(face_index).BCProperty());
		  }

		  // User Information of the results of the operation
		  Vec3d ref_colour = new Vec3d(0.0, 1.0, 0.0);
		  PrintMessage(3, "Colour based Boundary Condition Property details:");
		  for (int i = 0; i < faces_sorted.Size(); i++)
		  {
			 if (colours_sorted[i] > 0)
			 {
				 ref_colour = all_colours.Elem(colours_sorted[i]);
			 }

			 PrintMessage(3, "BC Property: ", i + DefineConstants.DEFAULT_BCNUM);
			 PrintMessage(3, "   Nr. of Surface Elements = ", faces_sorted[i]);
			 PrintMessage(3, "   Colour Index = ", colours_sorted[i]);
			 PrintMessage(3, "   RGB Face Colour = ", ref_colour, "", "\n");
		  }
	   }

	  public static ostream operator << (ostream ost, MarkedTet mt)
	  {
		for (int i = 0; i < 4; i++)
		{
		  ost << mt.pnums[i] << " ";
		}

		ost << mt.matindex << " " << (int)mt.marked << " " << (int)mt.flagged << " " << (int)mt.tetedge1 << " " << (int)mt.tetedge2 << " ";

		ost << "faceedges = ";
		for (int i = 0; i < 4; i++)
		{
		  ost << (int)mt.faceedges[i] << " ";
		}

		ost << " order = ";
		ost << mt.incorder << " " << (int)mt.order << "\n";
		return ost;
	  }
	  public static istream operator >> (istream ost, MarkedTet mt)
	  {
		for (int i = 0; i < 4; i++)
		{
		  ost >> mt.pnums[i];
		}

		ost >> mt.matindex;

		int auxint;
		ost >> auxint;
		mt.marked = auxint;
		ost >> auxint;
		mt.flagged = auxint;
		ost >> auxint;
		mt.tetedge1 = auxint;
		ost >> auxint;
		mt.tetedge2 = auxint;

		char auxchar;

		for (int i = 0; i < 4; i++)
		{
		ost >> auxchar;
		mt.faceedges = StringFunctions.ChangeCharacter(mt.faceedges, i, auxchar);
		}

		ost >> mt.incorder;
		ost >> auxint;
		mt.order = auxint;
		return ost;
	  }


	  public static ostream operator << (ostream ost, MarkedPrism mp)
	  {
		for (int i = 0; i < 6; i++)
		{
		  ost << mp.pnums[i] << " ";
		}

		ost << mp.matindex << " " << mp.marked << " " << mp.markededge << " " << mp.incorder << " " << (int)mp.order << "\n";
		return ost;
	  }
	  public static istream operator >> (istream ist, MarkedPrism mp)
	  {
		for (int i = 0; i < 6; i++)
		{
		  ist >> mp.pnums[i];
		}

		ist >> mp.matindex >> mp.marked >> mp.markededge >> mp.incorder;
		int auxint;
		ist >> auxint;
		mp.order = auxint;
		return ist;
	  }


	  public static ostream operator << (ostream ost, MarkedIdentification mi)
	  {
		ost << mi.np << " ";
		for (int i = 0; i < 2 * mi.np; i++)
		{
		  ost << mi.pnums[i] << " ";
		}
		ost << mi.marked << " " << mi.markededge << " " << mi.incorder << " " << (int)mi.order << "\n";
		return ost;
	  }
	  public static istream operator >> (istream ist, MarkedIdentification mi)
	  {
		ist >> mi.np;
		for (int i = 0; i < 2 * mi.np; i++)
		{
		  ist >> mi.pnums[i];
		}
		ist >> mi.marked >> mi.markededge >> mi.incorder;
		int auxint;
		ist >> auxint;
		mi.order = auxint;
		return ist;
	  }

	  public static ostream operator << (ostream ost, MarkedTri mt)
	  {
		for (int i = 0; i < 3; i++)
		{
		  ost << mt.pnums[i] << " ";
		}
		for (int i = 0; i < 3; i++)
		{
		  ost << mt.pgeominfo[i] << " ";
		}
		ost << mt.marked << " " << mt.markededge << " " << mt.surfid << " " << mt.incorder << " " << (int)mt.order << "\n";
		return ost;
	  }
	  public static istream operator >> (istream ist, MarkedTri mt)
	  {
		for (int i = 0; i < 3; i++)
		{
		  ist >> mt.pnums[i];
		}
		for (int i = 0; i < 3; i++)
		{
		  ist >> mt.pgeominfo[i];
		}
		ist >> mt.marked >> mt.markededge >> mt.surfid >> mt.incorder;
		int auxint;
		ist >> auxint;
		mt.order = auxint;
		return ist;
	  }

	  public static ostream operator << (ostream ost, MarkedQuad mt)
	  {
		for (int i = 0; i < 4; i++)
		{
		  ost << mt.pnums[i] << " ";
		}
		for (int i = 0; i < 4; i++)
		{
		  ost << mt.pgeominfo[i] << " ";
		}
		ost << mt.marked << " " << mt.markededge << " " << mt.surfid << " " << mt.incorder << " " << (int)mt.order << "\n";
		return ost;
	  }
	  public static istream operator >> (istream ist, MarkedQuad mt)
	  {
		for (int i = 0; i < 4; i++)
		{
		  ist >> mt.pnums[i];
		}
		for (int i = 0; i < 4; i++)
		{
		  ist >> mt.pgeominfo[i];
		}
		ist >> mt.marked >> mt.markededge >> mt.surfid >> mt.incorder;
		int auxint;
		ist >> auxint;
		mt.order = auxint;
		return ist;
	  }




	  public static void PrettyPrint(ostream ost, MarkedTet mt)
	  {
		int te1 = mt.tetedge1;
		int te2 = mt.tetedge2;
		int order = mt.order;

		ost << "MT: " << mt.pnums[0] << " - " << mt.pnums[1] << " - " << mt.pnums[2] << " - " << mt.pnums[3] << "\n" << "marked edge: " << te1 << " - " << te2 << ", order = " << order << "\n";
		//for (int k = 0; k < 4; k++)
		//  ost << int(mt.faceedges[k]) << "  ";
		for (int k = 0; k < 4; k++)
		{
		ost << "face";
		for (int j = 0; j < 4; j++)
		{
		  if (j != k)
		  {
			ost << " " << mt.pnums[j];
		  }
		}
		for (int i = 0; i < 3; i++)
		{
		  for (int j = i + 1; j < 4; j++)
		  {
			if (i != k && j != k && (int)mt.faceedges[k] == 6 - k - i - j)
			{
			  ost << " marked edge " << mt.pnums[i] << " " << mt.pnums[j] << "\n";
			}
		  }
		}
		}
		ost << "\n";
	  }




	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int[][] BTSortEdges_tetedges =
	  {
		  new int[] {1, 2},
		  new int[] {1, 3},
		  new int[] {1, 4},
		  new int[] {2, 3},
		  new int[] {2, 4},
		  new int[] {3, 4}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int[][] BTSortEdges_prismedges =
	  {
		  new int[] {1, 2},
		  new int[] {1, 3},
		  new int[] {2, 3},
		  new int[] {4, 5},
		  new int[] {4, 6},
		  new int[] {5, 6},
		  new int[] {1, 4},
		  new int[] {2, 5},
		  new int[] {3, 6}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int[][] BTSortEdges_trigedges =
	  {
		  new int[] {1, 2},
		  new int[] {2, 3},
		  new int[] {3, 1}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int[][] BTSortEdges_quadedges =
	  {
		  new int[] {1, 2},
		  new int[] {2, 3},
		  new int[] {3, 4},
		  new int[] {4, 1}
	  };

	  public static int BTSortEdges(Mesh mesh, Array< Array<int,PointIndex.BASE> > idmaps, INDEX_2_CLOSED_HASHTABLE<int> edgenumber)
	  {
		PrintMessage(4, "sorting ... ");

		//  if (mesh.PureTetMesh())
		if (true)
		{
		// new, fast version

		Array<INDEX_2> edges = new Array<INDEX_2>();
		Array<int> eclasses = new Array<int>();

		int i;
		int j;
		int k;
		int cntedges = 0;
		int go_on;
		int ned = 0;

		// enumerate edges:
		for (i = 1; i <= mesh.GetNE(); i++)
		{
			Element el = mesh.VolumeElement(i);
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//		static int tetedges[6][2] = { { 1, 2 }, { 1, 3 }, { 1, 4 }, { 2, 3 }, { 2, 4 }, { 3, 4 } };
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//		static int prismedges[9][2] = { { 1, 2 }, { 1, 3 }, { 2, 3 }, { 4, 5 }, { 4, 6 }, { 5, 6 }, { 1, 4 }, { 2, 5 }, { 3, 6 } };
			int[][] pyramidedges =
			{
				new int[] {1, 2},
				new int[] {3, 4},
				new int[] {1, 5},
				new int[] {2, 5},
				new int[] {3, 5},
				new int[] {4, 5}
			};

			int[] tip = 0;

			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TET:
			  case ELEMENT_TYPE.TET10:
			  {
			  tip = BTSortEdges_tetedges;
			  ned = 6;
			  break;
			  }
			  case ELEMENT_TYPE.PRISM:
			  case ELEMENT_TYPE.PRISM12:
			  {
			  tip = BTSortEdges_prismedges;
			  ned = 6;
			  break;
			  }
			  case ELEMENT_TYPE.PYRAMID:
			  {
			  tip = pyramidedges;
			  ned = 6;
			  break;
			  }
				  default:
					throw new Exception("Bisect, element type not handled in switch");
			}

			for (j = 0; j < ned; j++)
			{
			INDEX_2 i2 = new INDEX_2(el.PNum(tip[j][0]), el.PNum(tip[j][1]));
			i2.Sort();
			//(*testout) << "edge " << i2 << endl;
			if (!edgenumber.Used(i2))
			{
				cntedges++;
				edges.Append(i2);
				edgenumber.Set(i2, cntedges);
			}
			}
		}

		// additional surface edges:
		for (i = 1; i <= mesh.GetNSE(); i++)
		{
			Element2d el = mesh.SurfaceElement(i);
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//		static int trigedges[3][2] = { { 1, 2 }, { 2, 3 }, { 3, 1 } };

	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//		static int quadedges[4][2] = { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 1 } };


			int[] tip = 0;

			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TRIG:
			  case ELEMENT_TYPE.TRIG6:
			  {
			  tip = BTSortEdges_trigedges;
			  ned = 3;
			  break;
			  }
			  case ELEMENT_TYPE.QUAD:
			  case ELEMENT_TYPE.QUAD6:
			  {
			  tip = BTSortEdges_quadedges;
			  ned = 4;
			  break;
			  }
			  default:
			  {
			  cerr << "Error: Sort for Bisect, SE has " << el.GetNP() << " points" << "\n";
			  ned = 0;
			  }
		break;
			}

			for (j = 0; j < ned; j++)
			{
			INDEX_2 i2 = new INDEX_2(el.PNum(tip[j][0]), el.PNum(tip[j][1]));
			i2.Sort();
			if (!edgenumber.Used(i2))
			{
				cntedges++;
				edges.Append(i2);
				edgenumber.Set(i2, cntedges);
			}
			}
		}





		eclasses.SetSize(cntedges);
		for (i = 1; i <= cntedges; i++)
		{
		  eclasses.Elem(i) = i;
		}

		// identify edges in element stack
		do
		{
			go_on = 0;
			for (i = 1; i <= mesh.GetNE(); i++)
			{
			Element el = mesh.VolumeElement(i);

			if (el.GetType() != ELEMENT_TYPE.PRISM && el.GetType() != ELEMENT_TYPE.PRISM12 && el.GetType() != ELEMENT_TYPE.PYRAMID)
			{
			  continue;
			}

			int[][] prismpairs =
			{
				new int[] {1, 2, 4, 5},
				new int[] {2, 3, 5, 6},
				new int[] {1, 3, 4, 6}
			};

			int[][] pyramidpairs =
			{
				new int[] {1, 2, 4, 3},
				new int[] {1, 5, 4, 5},
				new int[] {2, 5, 3, 5}
			};

			int[] pairs = 0;
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.PRISM:
			  case ELEMENT_TYPE.PRISM12:
			  {
				  pairs = prismpairs;
				  break;
			  }
			  case ELEMENT_TYPE.PYRAMID:
			  {
				  pairs = pyramidpairs;
				  break;
			  }
					  default:
						throw new Exception("Bisect, element type not handled in switch, 2");
			}

			for (j = 0; j < 3; j++)
			{
				INDEX_2 e1 = new INDEX_2(el.PNum(pairs[j][0]), el.PNum(pairs[j][1]));
				INDEX_2 e2 = new INDEX_2(el.PNum(pairs[j][2]), el.PNum(pairs[j][3]));
				e1.Sort();
				e2.Sort();

				int eclass1 = edgenumber.Get(e1);
				int eclass2 = edgenumber.Get(e2);

				//		  (*testout) << "identify edges " << eclass1 << "-" << eclass2 << endl;

				if (eclasses.Get(eclass1) > eclasses.Get(eclass2))
				{
				eclasses.Elem(eclass1) = eclasses.Get(eclass2);
				go_on = 1;
				}
				else if (eclasses.Get(eclass2) > eclasses.Get(eclass1))
				{
				eclasses.Elem(eclass2) = eclasses.Get(eclass1);
				go_on = 1;
				}
			}
			}

			for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
			{
			Element2d el2d = mesh[sei];

			for (i = 0; i < el2d.GetNP(); i++)
			{
				INDEX_2 e1 = new INDEX_2(el2d[i], el2d[(i + 1) % el2d.GetNP()]);
				e1.Sort();
				INDEX_2 e2 = new INDEX_2();

				for (k = 0; k < idmaps.Size(); k++)
				{
				e2.I1() = (*idmaps[k])[e1.I1()];
				e2.I2() = (*idmaps[k])[e1.I2()];

				if (e2.I1() == 0 || e2.I2() == 0 || e1.I1() == e2.I1() || e1.I2() == e2.I2())
				{
				  continue;
				}

				e2.Sort();
				if (!edgenumber.Used(e2))
				{
				  continue;
				}


				int eclass1 = edgenumber.Get(e1);
				int eclass2 = edgenumber.Get(e2);

				if (eclasses.Get(eclass1) > eclasses.Get(eclass2))
				{
					eclasses.Elem(eclass1) = eclasses.Get(eclass2);


					go_on = 1;
				}
				else if (eclasses.Get(eclass2) > eclasses.Get(eclass1))
				{
					eclasses.Elem(eclass2) = eclasses.Get(eclass1);
					go_on = 1;
				}
				}
			}

			}

		} while (go_on != 0);

	// 	for (i = 1; i <= cntedges; i++)
	// 	  {
	// 	    (*testout) << "edge " << i << ": " 
	// 		       << edges.Get(i).I1() << "-" << edges.Get(i).I2()
	// 		       << ", class = " << eclasses.Get(i) << endl;
	// 	  }

		// compute classlength:
		Array<double> edgelength = new Array<double>(cntedges);

		/*
		for (i = 1; i <= cntedges; i++)
		  edgelength.Elem(i) = 1e20;
		*/

		for (i = 1; i <= cntedges; i++)
		{
			INDEX_2 edge = edges.Get(i);
			double elen = Dist(new mesh.Point(edge.I1()), new mesh.Point(edge.I2()));
			edgelength.Elem(i) = elen;
		}

		/*
		  for (i = 1; i <= mesh.GetNE(); i++)
		  {
		  const Element & el = mesh.VolumeElement (i);
		  
		  if (el.GetType() == TET)
		  {
		  for (j = 1; j <= 3; j++)
		  for (k = j+1; k <= 4; k++)
		  {
		  INDEX_2 i2(el.PNum(j), el.PNum(k));
		  i2.Sort();
			    
		  int enr = edgenumber.Get(i2);
		  double elen = Dist (mesh.Point (i2.I1()), mesh.Point (i2.I2()));
		  if (elen < edgelength.Get(enr))
		  edgelength.Set (enr, elen);
		  }
		  }
		  else if (el.GetType() == PRISM)
		  {
		  for (j = 1; j <= 3; j++)
		  {
		  k = (j % 3) + 1;
			  
		  INDEX_2 i2(el.PNum(j), el.PNum(k));
		  i2.Sort();
			  
		  int enr = edgenumber.Get(i2);
		  double elen = Dist (mesh.Point (i2.I1()), mesh.Point (i2.I2()));
		  if (elen < edgelength.Get(enr))
		  edgelength.Set (enr, elen);
			  
		  i2 = INDEX_2(el.PNum(j+3), el.PNum(k+3));
		  i2.Sort();
			  
		  enr = edgenumber.Get(i2);
		  elen = Dist (mesh.Point (i2.I1()), mesh.Point (i2.I2()));
		  if (elen < edgelength.Get(enr))
		  edgelength.Set (enr, elen);
			  
		  if (!edgenumber.Used(i2))
		  {
		  cntedges++;
		  edgenumber.Set(i2, cntedges);
		  }
		  i2 = INDEX_2(el.PNum(j), el.PNum(j+3));
		  i2.Sort();
			  
		  enr = edgenumber.Get(i2);
		  elen = Dist (mesh.Point (i2.I1()), mesh.Point (i2.I2()));
		  if (elen < edgelength.Get(enr))
		  edgelength.Set (enr, elen);
		  }
		  }
		  }
		*/


		for (i = 1; i <= cntedges; i++)
		{
			if (eclasses.Get(i) != i)
			{
			if (edgelength.Get(i) < edgelength.Get(eclasses.Get(i)))
			{
			  edgelength.Elem(eclasses.Get(i)) = edgelength.Get(i);
			}
			edgelength.Elem(i) = 1e20;
			}
		}


		TABLE<int> eclasstab = new TABLE<int>(cntedges);
		for (i = 1; i <= cntedges; i++)
		{
		  eclasstab.Add1(eclasses.Get(i), i);
		}


		// sort edges:
		Array<int> sorted = new Array<int>(cntedges);

		QuickSort(edgelength, sorted);

		int cnt = 0;
		for (i = 1; i <= cntedges; i++)
		{
			int ii = sorted.Get(i);
			for (j = 1; j <= eclasstab.EntrySize(ii); j++)
			{
			cnt++;
			edgenumber.Set(edges.Get(eclasstab.Get(ii, j)), cnt);
			}
		}
		return cnt;
		}

		else

		{
		// old version

		int i;
		int j;
		int cnt = 0;
		int found;
		double len2;
		double maxlen2;
		INDEX_2 ep = new INDEX_2();

		// sort edges by length, parallel edges (on prisms)
		// are added in blocks

		do
		{
			found = 0;
			maxlen2 = 1e30;

			for (i = 1; i <= mesh.GetNE(); i++)
			{
			Element el = mesh.VolumeElement(i);
			int ned;
			int[][] BTSortEdges_tetedges =
			{
				new int[] {1, 2},
				new int[] {1, 3},
				new int[] {1, 4},
				new int[] {2, 3},
				new int[] {2, 4},
				new int[] {3, 4}
			};
			int[][] BTSortEdges_prismedges =
			{
				new int[] {1, 2},
				new int[] {1, 3},
				new int[] {2, 4},
				new int[] {4, 5},
				new int[] {4, 6},
				new int[] {5, 6}
			};
			int[][] pyramidedges =
			{
				new int[] {1, 2},
				new int[] {3, 4},
				new int[] {1, 5},
				new int[] {2, 5},
				new int[] {3, 5},
				new int[] {4, 5}
			};

			int[] tip = new int[2];

			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TET:
			  {
				  tip = BTSortEdges_tetedges;
				  ned = 6;
				  break;
			  }
			  case ELEMENT_TYPE.PRISM:
			  {
				  tip = BTSortEdges_prismedges;
				  ned = 6;
				  break;
			  }
			  case ELEMENT_TYPE.PYRAMID:
			  {
				  tip = pyramidedges;
				  ned = 6;
				  break;
			  }
					  default:
						throw new Exception("Bisect, element type not handled in switch, 3");
			}

			for (j = 0; j < ned; j++)
			{
				INDEX_2 i2 = new INDEX_2(el.PNum(tip[j][0]), el.PNum(tip[j][1]));
				i2.Sort();
				if (!edgenumber.Used(i2))
				{
				len2 = Dist(new mesh.Point(i2.I1()), new mesh.Point(i2.I2()));
				if (len2 < maxlen2)
				{
					maxlen2 = len2;
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: ep = i2;
					ep.CopyFrom(i2);
					found = 1;
				}
				}
			}
			}
			if (found != 0)
			{
			cnt++;
			edgenumber.Set(ep, cnt);


			// find connected edges:
			int go_on = 0;
			do
			{
				go_on = 0;
				for (i = 1; i <= mesh.GetNE(); i++)
				{
				Element el = mesh.VolumeElement(i);
				if (el.GetNP() != 6)
				{
					continue;
				}

				int[][] prismpairs =
				{
					new int[] {1, 2, 4, 5},
					new int[] {2, 3, 5, 6},
					new int[] {1, 3, 4, 6}
				};

				int[][] pyramidpairs =
				{
					new int[] {1, 2, 4, 3},
					new int[] {1, 5, 4, 5},
					new int[] {2, 5, 3, 5}
				};

				int[] pairs = new int[4];
				switch (el.GetType())
				{
				  case ELEMENT_TYPE.PRISM:
				  {
					  pairs = prismpairs;
					  break;
				  }
				  case ELEMENT_TYPE.PYRAMID:
				  {
					  pairs = pyramidpairs;
					  break;
				  }
							  default:
								throw new Exception("Bisect, element type not handled in switch, 3a");
				}

				for (j = 0; j < 3; j++)
				{
					INDEX_2 e1 = new INDEX_2(el.PNum(pairs[j][0]), el.PNum(pairs[j][1]));
					INDEX_2 e2 = new INDEX_2(el.PNum(pairs[j][2]), el.PNum(pairs[j][3]));
					e1.Sort();
					e2.Sort();

					int used1 = edgenumber.Used(e1);
					int used2 = edgenumber.Used(e2);

					if (used1 != 0 && used2 == 0)
					{
					cnt++;
					edgenumber.Set(e2, cnt);
					go_on = 1;
					}
					if (used2 != 0 && used1 == 0)
					{
					cnt++;
					edgenumber.Set(e1, cnt);
					go_on = 1;
					}
				}
				}
			} while (go_on != 0);
			}
		} while (found != 0);

		return cnt;
		}
	  }




	  public static void BTDefineMarkedTet(Element el, INDEX_2_CLOSED_HASHTABLE<int> edgenumber, MarkedTet mt)
	  {
		for (int i = 0; i < 4; i++)
		{
		  mt.pnums[i] = el[i];
		}

		mt.marked = 0;
		mt.flagged = 0;

		mt.incorder = false;
		mt.order = 1;

		int val = 0;
		// find marked edge of tet:
		for (int i = 0; i < 3; i++)
		{
		  for (int j = i + 1; j < 4; j++)
		  {
		  INDEX_2 i2 = new INDEX_2(mt.pnums[i], mt.pnums[j]);
		  i2.Sort();
		  int hval = edgenumber.Get(i2);
		  if (hval > val)
		  {
			  val = hval;
			  mt.tetedge1 = i;
			  mt.tetedge2 = j;
		  }
		  }
		}


		// find marked edges of faces:
		for (int k = 0; k < 4; k++)
		{
		val = 0;
		for (int i = 0; i < 3; i++)
		{
		  for (int j = i + 1; j < 4; j++)
		  {
			if (i != k && j != k)
			{
			INDEX_2 i2 = new INDEX_2(mt.pnums[i], mt.pnums[j]);
			i2.Sort();
			int hval = edgenumber.Get(i2);
			if (hval > val)
			{
				val = hval;
						int hi = 6 - k - i - j;
						mt.faceedges = StringFunctions.ChangeCharacter(mt.faceedges, k, (char)hi);
			}
			}
		  }
		}
		}
	  }




	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int[] BTDefineMarkedPrism_map = {1, 2, 5, 4, 3, 5};
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int[] BTDefineMarkedPrism_map = {1, 4, 3, 2, 4, 3};

	  public static void BTDefineMarkedPrism(Element el, INDEX_2_CLOSED_HASHTABLE<int> edgenumber, MarkedPrism mp)
	  {
		if (el.GetType() == ELEMENT_TYPE.PRISM || el.GetType() == ELEMENT_TYPE.PRISM12)
		{
		for (int i = 0; i < 6; i++)
		{
		  mp.pnums[i] = el[i];
		}
		}
		else if (el.GetType() == ELEMENT_TYPE.PYRAMID)
		{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static int map[6] = { 1, 2, 5, 4, 3, 5 };
		for (int i = 0; i < 6; i++)
		{
		  mp.pnums[i] = el.PNum(BTDefineMarkedPrism_map[i]);
		}
		}
		else if (el.GetType() == ELEMENT_TYPE.TET || el.GetType() == ELEMENT_TYPE.TET10)
		{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static int map[6] = { 1, 4, 3, 2, 4, 3 };
		for (int i = 0; i < 6; i++)
		{
		  mp.pnums[i] = el.PNum(BTDefineMarkedPrism_map[i]);
		}

		}
		else
		{
		PrintSysError("Define marked prism called for non-prism and non-pyramid");
		}



		mp.marked = 0;

		mp.incorder = false;
		mp.order = 1;

		int val = 0;
		for (int i = 0; i < 2; i++)
		{
		  for (int j = i + 1; j < 3; j++)
		  {
		  INDEX_2 i2 = new INDEX_2(mp.pnums[i], mp.pnums[j]);
		  i2.Sort();
		  int hval = edgenumber.Get(i2);
		  if (hval > val)
		  {
			  val = hval;
			  mp.markededge = 3 - i - j;
		  }
		  }
		}
	  }



	  public static bool BTDefineMarkedId(Element2d el, INDEX_2_CLOSED_HASHTABLE<int> edgenumber, Array<int,PointIndex.BASE> idmap, MarkedIdentification mi)
	  {

		bool identified = true;
		mi.np = el.GetNP();
		int min1 = 0;
		int min2 = 0;
		for (int j = 0; identified && j < mi.np; j++)
		{
		mi.pnums[j] = el[j];
		mi.pnums[j + mi.np] = idmap[el[j]];

		if (j == 0 || el[j] < min1)
		{
		  min1 = el[j];
		}
		if (j == 0 || mi.pnums[j + mi.np] < min2)
		{
		  min2 = mi.pnums[j + mi.np];
		}

		identified = (mi.pnums[j + mi.np] != 0 && mi.pnums[j + mi.np] != mi.pnums[j]);
		}

		identified = identified && (min1 < min2);

		if (identified)
		{
		mi.marked = 0;

		mi.incorder = false;
		mi.order = 1;

		int val = 0;
		for (int i = 0; i < mi.np; i++)
		{
			INDEX_2 i2 = new INDEX_2(mi.pnums[i], mi.pnums[(i + 1) % mi.np]);
			i2.Sort();
			int hval = edgenumber.Get(i2);
			if (hval > val)
			{
			val = hval;
			mi.markededge = i;
			}
		}
		}

		return identified;
	  }


	  public static void BTDefineMarkedTri(Element2d el, INDEX_2_CLOSED_HASHTABLE<int> edgenumber, MarkedTri mt)
	  {
		for (int i = 0; i < 3; i++)
		{
		mt.pnums[i] = el[i];
		mt.pgeominfo[i] = el.GeomInfoPi(i + 1);
		}

		mt.marked = 0;
		mt.surfid = el.GetIndex();

		mt.incorder = false;
		mt.order = 1;

		int val = 0;
		for (int i = 0; i < 2; i++)
		{
		  for (int j = i + 1; j < 3; j++)
		  {
		  INDEX_2 i2 = new INDEX_2(mt.pnums[i], mt.pnums[j]);
		  i2.Sort();
		  int hval = edgenumber.Get(i2);
		  if (hval > val)
		  {
			  val = hval;
			  mt.markededge = 3 - i - j;
		  }
		  }
		}
	  }



	  public static void PrettyPrint(ostream ost, MarkedTri mt)
	  {
		ost << "MarkedTrig: " << "\n";
		ost << "  pnums = ";
		for (int i = 0; i < 3; i++)
		{
			ost << mt.pnums[i] << " ";
		}
		ost << "\n";
		ost << "  marked = " << mt.marked << ", markededge=" << mt.markededge << "\n";
		for (int i = 0; i < 2; i++)
		{
		  for (int j = i + 1; j < 3; j++)
		  {
		if (mt.markededge == 3 - i - j)
		{
		  ost << "  marked edge pnums = " << mt.pnums[i] << " " << mt.pnums[j] << "\n";
		}
		  }
		}
	  }


	  public static void PrettyPrint(ostream ost, MarkedQuad mq)
	  {
		ost << "MarkedQuad: " << "\n";
		ost << "  pnums = ";
		for (int i = 0; i < 4; i++)
		{
			ost << mq.pnums[i] << " ";
		}
		ost << "\n";
		ost << "  marked = " << mq.marked << ", markededge=" << mq.markededge << "\n";
	  }





	  public static void BTDefineMarkedQuad(Element2d el, INDEX_2_CLOSED_HASHTABLE<int> edgenumber, MarkedQuad mq)
	  {
		for (int i = 0; i < 4; i++)
		{
		  mq.pnums[i] = el[i];
		}
		Swap(ref mq.pnums[2], ref mq.pnums[3]);

		mq.marked = 0;
		mq.markededge = 0;
		mq.surfid = el.GetIndex();
	  }




	  // mark elements due to local h
	  public static int BTMarkTets(Array<MarkedTet> mtets, Array<MarkedPrism> mprisms, Mesh mesh)
	  {
		int marked = 0;

		int np = mesh.GetNP();
		Vector hv = new Vector(np);
		for (int i = 0; i < np; i++)
		{
		  hv(i) = mesh.GetH(new mesh.Point(i + 1));
		}

		double hfac = 1;

		for (int step = 1; step <= 2; step++)
		{
		for (int i = 1; i <= mtets.Size(); i++)
		{
			double h = 0;

			for (int j = 0; j < 3; j++)
			{
			  for (int k = j + 1; k < 4; k++)
			  {
			  const Point < 3> & p1 = new mesh.Point(mtets.Get(i).pnums[j]);
			  const Point < 3> & p2 = new mesh.Point(mtets.Get(i).pnums[k]);
			  double hh = Dist2(p1, p2);
			  if (hh > h)
			  {
				  h = hh;
			  }
			  }
			}
			h = ngsimd.GlobalMembers.sqrt(h);

			double hshould = 1e10;
			for (int j = 0; j < 4; j++)
			{
			double hi = hv(mtets.Get(i).pnums[j] - 1);
			if (hi < hshould)
			{
			  hshould = hi;
			}
			}


			if (step == 1)
			{
			if (h / hshould > hfac)
			{
			  hfac = h / hshould;
			}
			}
			else
			{
			if (h > hshould * hfac)
			{
				mtets.Elem(i).marked = 1;
				marked = 1;
			}
			else
			{
			  mtets.Elem(i).marked = 0;
			}
			}

		}
		for (int i = 1; i <= mprisms.Size(); i++)
		{
			double h = 0;

			for (int j = 0; j < 2; j++)
			{
			  for (int k = j + 1; k < 3; k++)
			  {
			  const Point < 3> & p1 = new mesh.Point(mprisms.Get(i).pnums[j]);
			  const Point < 3> & p2 = new mesh.Point(mprisms.Get(i).pnums[k]);
			  double hh = Dist2(p1, p2);
			  if (hh > h)
			  {
				  h = hh;
			  }
			  }
			}
			h = ngsimd.GlobalMembers.sqrt(h);

			double hshould = 1e10;
			for (int j = 0; j < 6; j++)
			{
			double hi = hv(mprisms.Get(i).pnums[j] - 1);
			if (hi < hshould)
			{
			  hshould = hi;
			}
			}


			if (step == 1)
			{
			if (h / hshould > hfac)
			{
			  hfac = h / hshould;
			}
			}
			else
			{
			if (h > hshould * hfac)
			{
				mprisms.Elem(i).marked = 1;
				marked = 1;
			}
			else
			{
			  mprisms.Elem(i).marked = 0;
			}
			}

		}



		if (step == 1)
		{
			if (hfac > 2)
			{
			  hfac /= 2;
			}
			else
			{
			  hfac = 1;
			}
		}

		}
		return marked;
	  }














	  public static void BTBisectTet(MarkedTet oldtet, int newp, MarkedTet newtet1, MarkedTet newtet2)
	  {
	#if DEBUG
		*testout << "bisect tet " << oldtet << "\n";
	#endif


		// points vis a vis from tet-edge
		int vis1;
		int vis2;
		vis1 = 0;
		while (vis1 == oldtet.tetedge1 || vis1 == oldtet.tetedge2)
		{
		  vis1++;
		}
		vis2 = 6 - vis1 - oldtet.tetedge1 - oldtet.tetedge2;





		// is tet of type P ?
		int istypep = 0;
		for (int i = 0; i < 4; i++)
		{
		int cnt = 0;
		for (int j = 0; j < 4; j++)
		{
		  if (oldtet.faceedges[j] == i)
		  {
			cnt++;
		  }
		}
		if (cnt == 3)
		{
		  istypep = 1;
		}
		}



		for (int i = 0; i < 4; i++)
		{
		newtet1.pnums[i] = oldtet.pnums[i];
		newtet2.pnums[i] = oldtet.pnums[i];
		}
		newtet1.flagged = istypep && !oldtet.flagged;
		newtet2.flagged = istypep && !oldtet.flagged;

		int nm = oldtet.marked - 1;
		if (nm < 0)
		{
			nm = 0;
		}
		newtet1.marked = nm;
		newtet2.marked = nm;

	#if DEBUG
		*testout << "newtet1,before = " << newtet1 << "\n";
		*testout << "newtet2,before = " << newtet2 << "\n";
	#endif

		for (int i = 0; i < 4; i++)
		{
		if (i == oldtet.tetedge1)
		{
			newtet2.pnums[i] = newp;
			newtet2.faceedges = StringFunctions.ChangeCharacter(newtet2.faceedges, i, oldtet.faceedges[i]); // inherited face
			newtet2.faceedges = StringFunctions.ChangeCharacter(newtet2.faceedges, vis1, i); // cut faces
			newtet2.faceedges = StringFunctions.ChangeCharacter(newtet2.faceedges, vis2, i);

			int j = 0;
			while (j == i || j == oldtet.faceedges[i])
			{
			  j++;
			}
			int k = 6 - i - oldtet.faceedges[i] - j;
			newtet2.tetedge1 = j; // tet-edge
			newtet2.tetedge2 = k;

			// new face:
			if (istypep != 0 && oldtet.flagged)
			{
					int hi = 6 - oldtet.tetedge1 - j - k;
					newtet2.faceedges = StringFunctions.ChangeCharacter(newtet2.faceedges, oldtet.tetedge2, (char)hi);
			}
			else
			{
			  newtet2.faceedges = StringFunctions.ChangeCharacter(newtet2.faceedges, oldtet.tetedge2, oldtet.tetedge1);
			}

	#if DEBUG
				*testout << "i = " << i << ", j = " << j << " k = " << k << " oldtet.tetedge1 = " << oldtet.tetedge1 << " oldtet.tetedge2 = " << oldtet.tetedge2 << "   6-oldtet.tetedge1-j-k = " << 6 - oldtet.tetedge1 - j - k << "   6-oldtet.tetedge1-j-k = " << (short)(6 - oldtet.tetedge1 - j - k) << "\n";
				*testout << "vis1 = " << vis1 << ", vis2 = " << vis2 << "\n";
				for (int j = 0; j < 4; j++)
				{
				  if (newtet2.faceedges[j] > 3)
				  {
					  *testout << "ERROR1" << "\n";
				  }
				}
	#endif
		}

		if (i == oldtet.tetedge2)
		{
			newtet1.pnums[i] = newp;
			newtet1.faceedges = StringFunctions.ChangeCharacter(newtet1.faceedges, i, oldtet.faceedges[i]); // inherited face
			newtet1.faceedges = StringFunctions.ChangeCharacter(newtet1.faceedges, vis1, i);
			newtet1.faceedges = StringFunctions.ChangeCharacter(newtet1.faceedges, vis2, i);
			int j = 0;
			while (j == i || j == oldtet.faceedges[i])
			{
			  j++;
			}
			int k = 6 - i - oldtet.faceedges[i] - j;
			newtet1.tetedge1 = j;
			newtet1.tetedge2 = k;

			// new face:
			if (istypep != 0 && oldtet.flagged)
			{
					int hi = 6 - oldtet.tetedge2 - j - k;
					newtet1.faceedges = StringFunctions.ChangeCharacter(newtet1.faceedges, oldtet.tetedge1, (char)hi);
			}
			else
			{
			  newtet1.faceedges = StringFunctions.ChangeCharacter(newtet1.faceedges, oldtet.tetedge1, oldtet.tetedge2);
			}

	#if DEBUG
				for (int j = 0; j < 4; j++)
				{
				  if (newtet2.faceedges[j] > 3)
				  {
					  *testout << "ERROR2" << "\n";
				  }
				}
	#endif
		}
		}

		newtet1.matindex = oldtet.matindex;
		newtet2.matindex = oldtet.matindex;
		newtet1.incorder = false;
		newtet1.order = oldtet.order;
		newtet2.incorder = false;
		newtet2.order = oldtet.order;

		// *testout << "newtet1 =  " << newtet1 << endl;
		// *testout << "newtet2 =  " << newtet2 << endl;
	  }




	  public static void BTBisectPrism(MarkedPrism oldprism, int newp1, int newp2, MarkedPrism newprism1, MarkedPrism newprism2)
	  {
		for (int i = 0; i < 6; i++)
		{
		newprism1.pnums[i] = oldprism.pnums[i];
		newprism2.pnums[i] = oldprism.pnums[i];
		}

		int pe1 = 0;
		if (pe1 == oldprism.markededge)
		{
		  pe1++;
		}
		int pe2 = 3 - oldprism.markededge - pe1;

		newprism1.pnums[pe2] = newp1;
		newprism1.pnums[pe2 + 3] = newp2;
		newprism1.markededge = pe2;
		newprism2.pnums[pe1] = newp1;
		newprism2.pnums[pe1 + 3] = newp2;
		newprism2.markededge = pe1;

		newprism1.matindex = oldprism.matindex;
		newprism2.matindex = oldprism.matindex;

		int nm = oldprism.marked - 1;
		if (nm < 0)
		{
			nm = 0;
		}
		newprism1.marked = nm;
		newprism2.marked = nm;

		newprism1.incorder = false;
		newprism1.order = oldprism.order;
		newprism2.incorder = false;
		newprism2.order = oldprism.order;
	  }


	  public static void BTBisectIdentification(MarkedIdentification oldid, Array<PointIndex> newp, MarkedIdentification newid1, MarkedIdentification newid2)
	  {
		for (int i = 0; i < 2 * oldid.np; i++)
		{
		newid1.pnums[i] = oldid.pnums[i];
		newid2.pnums[i] = oldid.pnums[i];
		}
		newid1.np = newid2.np = oldid.np;

		if (oldid.np == 2)
		{
			newid1.pnums[1] = newp[0];
			newid2.pnums[0] = newp[0];
			newid1.pnums[3] = newp[1];
			newid2.pnums[2] = newp[1];
			newid1.markededge = 0;
			newid2.markededge = 0;
		}

		if (oldid.np == 3)
		{
		newid1.pnums[(oldid.markededge+1) % 3] = newp[0];
		newid1.pnums[(oldid.markededge+1) % 3 + 3] = newp[1];
		newid1.markededge = (oldid.markededge+2) % 3;

		newid2.pnums[oldid.markededge] = newp[0];
		newid2.pnums[oldid.markededge+3] = newp[1];
		newid2.markededge = (oldid.markededge+1) % 3;
		}
		else if (oldid.np == 4)
		{
		newid1.pnums[(oldid.markededge+1) % 4] = newp[0];
		newid1.pnums[(oldid.markededge+2) % 4] = newp[2];
		newid1.pnums[(oldid.markededge+1) % 4 + 4] = newp[1];
		newid1.pnums[(oldid.markededge+2) % 4 + 4] = newp[3];
		newid1.markededge = (oldid.markededge+3) % 4;

		newid2.pnums[oldid.markededge] = newp[0];
		newid2.pnums[(oldid.markededge+3) % 4] = newp[2];
		newid2.pnums[oldid.markededge+4] = newp[1];
		newid2.pnums[(oldid.markededge+3) % 4 + 4] = newp[3];
		newid2.markededge = (oldid.markededge+1) % 4;
		}


		int nm = oldid.marked - 1;
		if (nm < 0)
		{
			nm = 0;
		}
		newid1.marked = newid2.marked = nm;

		newid1.incorder = newid2.incorder = false;
		newid1.order = newid2.order = oldid.order;
	  }



	  public static void BTBisectTri(MarkedTri oldtri, int newp, PointGeomInfo newpgi, MarkedTri newtri1, MarkedTri newtri2)
	  {
		for (int i = 0; i < 3; i++)
		{
		newtri1.pnums[i] = oldtri.pnums[i];
		newtri1.pgeominfo[i] = oldtri.pgeominfo[i];
		newtri2.pnums[i] = oldtri.pnums[i];
		newtri2.pgeominfo[i] = oldtri.pgeominfo[i];
		}

		int pe1 = 0;
		if (pe1 == oldtri.markededge)
		{
		  pe1++;
		}
		int pe2 = 3 - oldtri.markededge - pe1;

		newtri1.pnums[pe2] = newp;
		newtri1.pgeominfo[pe2] = newpgi;
		newtri1.markededge = pe2;

		newtri2.pnums[pe1] = newp;
		newtri2.pgeominfo[pe1] = newpgi;
		newtri2.markededge = pe1;


		newtri1.surfid = oldtri.surfid;
		newtri2.surfid = oldtri.surfid;

		int nm = oldtri.marked - 1;
		if (nm < 0)
		{
			nm = 0;
		}
		newtri1.marked = nm;
		newtri2.marked = nm;

		newtri1.incorder = false;
		newtri1.order = oldtri.order;
		newtri2.incorder = false;
		newtri2.order = oldtri.order;
	  }


	  public static void BTBisectQuad(MarkedQuad oldquad, int newp1, PointGeomInfo npgi1, int newp2, PointGeomInfo npgi2, MarkedQuad newquad1, MarkedQuad newquad2)
	  {
		for (int i = 0; i < 4; i++)
		{
		newquad1.pnums[i] = oldquad.pnums[i];
		newquad1.pgeominfo[i] = oldquad.pgeominfo[i];
		newquad2.pnums[i] = oldquad.pnums[i];
		newquad2.pgeominfo[i] = oldquad.pgeominfo[i];
		}

	/*    if (oldquad.marked==1) // he/sz: 2d quads or 3d prism
	    {   
	      newquad1.pnums[1] = newp1;
	      newquad1.pgeominfo[1] = npgi1;
	      newquad1.pnums[3] = newp2;
	      newquad1.pgeominfo[3] = npgi2;
	
	      newquad2.pnums[0] = newp1;
	      newquad2.pgeominfo[0] = npgi1;
	      newquad2.pnums[2] = newp2;
	      newquad2.pgeominfo[2] = npgi2;
	    }
	      
	    else if (oldquad.marked==2) // he/sz: 2d quads only
	    {
	      newquad1.pnums[0] = newp1;
	      newquad1.pnums[1] = newp2;
	      newquad1.pnums[3] = oldquad.pnums[2];  
	      newquad1.pnums[2] = oldquad.pnums[0]; 
	      newquad1.pgeominfo[0] = npgi1;
	      newquad1.pgeominfo[1] = npgi2;
	      newquad1.pgeominfo[3] = oldquad.pgeominfo[2]; 
	      newquad1.pgeominfo[2] = oldquad.pgeominfo[0];
	
	      newquad2.pnums[0] = newp2;
	      newquad2.pnums[1] = newp1;
	      newquad2.pnums[3] = oldquad.pnums[1];  
	      newquad2.pnums[2] = oldquad.pnums[3]; 
	      newquad2.pgeominfo[0] = npgi2;
	      newquad2.pgeominfo[1] = npgi1;
	      newquad2.pgeominfo[3] = oldquad.pgeominfo[1]; 
	      newquad2.pgeominfo[2] = oldquad.pgeominfo[3];
	    }
	      
	    */

		if (oldquad.markededge == 0 || oldquad.markededge == 2)
		{
		  newquad1.pnums[1] = newp1;
		  newquad1.pgeominfo[1] = npgi1;
		  newquad1.pnums[3] = newp2;
		  newquad1.pgeominfo[3] = npgi2;

		  newquad2.pnums[0] = newp1;
		  newquad2.pgeominfo[0] = npgi1;
		  newquad2.pnums[2] = newp2;
		  newquad2.pgeominfo[2] = npgi2;
		}
		else // 1 || 3
		{
		  newquad1.pnums[2] = newp1;
		  newquad1.pgeominfo[2] = npgi1;
		  newquad1.pnums[3] = newp2;
		  newquad1.pgeominfo[3] = npgi2;

		  newquad2.pnums[0] = newp1;
		  newquad2.pgeominfo[0] = npgi1;
		  newquad2.pnums[1] = newp2;
		  newquad2.pgeominfo[1] = npgi2;
		}
		newquad1.surfid = oldquad.surfid;
		newquad2.surfid = oldquad.surfid;

		int nm = oldquad.marked - 1;
		if (nm < 0)
		{
			nm = 0;
		}

		newquad1.marked = nm;
		newquad2.marked = nm;

		if (nm == 1)
		{
		  newquad1.markededge = 1;
		  newquad2.markededge = 1;
		}
		else
		{
		  newquad1.markededge = 0;
		  newquad2.markededge = 0;
		}

	  }


	  public static int MarkHangingIdentifications(Array<MarkedIdentification> mids, INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges)
	  {
		int hanging = 0;
		for (int i = 1; i <= mids.Size(); i++)
		{
		if (mids.Elem(i).marked)
		{
			hanging = 1;
			continue;
		}

		int np = mids.Get(i).np;
		for (int j = 0; j < np; j++)
		{
			INDEX_2 edge1 = new INDEX_2(mids.Get(i).pnums[j], mids.Get(i).pnums[(j + 1) % np]);
			INDEX_2 edge2 = new INDEX_2(mids.Get(i).pnums[j + np], mids.Get(i).pnums[((j + 1) % np) + np]);

			edge1.Sort();
			edge2.Sort();
			if (cutedges.Used(edge1) || cutedges.Used(edge2))
			{
			mids.Elem(i).marked = 1;
			hanging = 1;
			}
		}
		}

		return hanging;
	  }


	  /*
	  void IdentifyCutEdges(Mesh & mesh,
				INDEX_2_CLOSED_HASHTABLE<int> & cutedges)
	  {
	    int i,j,k;
	
	    Array< Array<int,PointIndex::BASE>* > idmaps;
	    for(i=1; i<=mesh.GetIdentifications().GetMaxNr(); i++)
	      {
		idmaps.Append(new Array<int,PointIndex::BASE>);
		mesh.GetIdentifications().GetMap(i,*idmaps.Last());
	      }
	
	
	    
	    for(SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	      {
		const Element2d & el2d = mesh[sei];
	
		for(i = 0; i < el2d.GetNP(); i++)
		  {
			INDEX_2 e1(el2d[i], el2d[(i+1) % el2d.GetNP()]);
			e1.Sort();
	
			if(!cutedges.Used(e1))
			  continue;
	
	
			for(k = 0; k < idmaps.Size(); k++)
			  {
			INDEX_2 e2((*idmaps[k])[e1.I1()],
				   (*idmaps[k])[e1.I2()]);
	
			if(e2.I1() == 0 || e2.I2() == 0 ||
			   e1.I1() == e2.I1() || e1.I2() == e2.I2())
			  continue;
	
			e2.Sort();
	
			if(cutedges.Used(e2))
			  continue;
	
			Point3d np = Center(mesh.Point(e2.I1()),
						mesh.Point(e2.I2()));
			int newp = mesh.AddPoint(np);
			cutedges.Set(e2,newp);
			(*testout) << "DAAA" << endl;
			  }
		  }
	      }
	
	    
	    for(i=0; i<idmaps.Size(); i++)
	      delete idmaps[i];
	    idmaps.DeleteAll();
	  }
	  */


	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int MarkHangingTets_timer = NgProfiler.CreateTimer("MarkHangingTets");

	  public static int MarkHangingTets(Array<MarkedTet> mtets, INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges, TaskManager tm)
	  {
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static int timer = NgProfiler::CreateTimer("MarkHangingTets");
		NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(MarkHangingTets_timer);

		int hanging = 0;
		// for (int i = 1; i <= mtets.Size(); i++)
	ParallelForRange(tm, mtets.Size(), (uint begin, uint end) =>
	{
		 bool my_hanging = false;
		 for (uint i = begin; i < end; i++)
		 {
			 MarkedTet teti = mtets[i];

			 if (teti.marked)
			 {
				 my_hanging = true;
				 continue;
			 }

			 for (int j = 0; j < 3; j++)
			 {
			   for (int k = j + 1; k < 4; k++)
			   {
				   INDEX_2 edge = new INDEX_2(teti.pnums[j], teti.pnums[k]);
				   edge.Sort();
				   if (cutedges.Used(edge))
				   {
					   teti.marked = 1;
					   my_hanging = true;
				   }
			   }
			 }
		 }
		 if (my_hanging)
		 {
			 hanging = true;
		 }
	});

		return hanging;
	  }



	  public static int MarkHangingPrisms(Array<MarkedPrism> mprisms, INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges)
	  {
		int hanging = 0;
		for (int i = 1; i <= mprisms.Size(); i++)
		{
		if (mprisms.Elem(i).marked)
		{
			hanging = 1;
			continue;
		}

		for (int j = 0; j < 2; j++)
		{
		  for (int k = j + 1; k < 3; k++)
		  {
			  INDEX_2 edge1 = new INDEX_2(mprisms.Get(i).pnums[j], mprisms.Get(i).pnums[k]);
			  INDEX_2 edge2 = new INDEX_2(mprisms.Get(i).pnums[j + 3], mprisms.Get(i).pnums[k + 3]);
			  edge1.Sort();
			  edge2.Sort();
			  if (cutedges.Used(edge1) || cutedges.Used(edge2))
			  {
			  mprisms.Elem(i).marked = 1;
			  hanging = 1;
			  }
		  }
		}
		}
		return hanging;
	  }



	  public static bool MarkHangingTris(Array<MarkedTri> mtris, INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges, TaskManager tm)
	  {
		bool hanging = false;
		// for (int i = 1; i <= mtris.Size(); i++)
		// for (auto & tri : mtris)
	ParallelForRange(tm, mtris.Size(), (uint begin, uint end) =>
	{
		 bool my_hanging = false;
		 for (uint i = begin; i < end; i++)
		 {
			 var tri = mtris[i];
			 if (tri.marked)
			 {
				 my_hanging = true;
				 continue;
			 }
			 for (int j = 0; j < 2; j++)
			 {
			   for (int k = j + 1; k < 3; k++)
			   {
				   INDEX_2 edge = new INDEX_2(tri.pnums[j], tri.pnums[k]);
				   edge.Sort();
				   if (cutedges.Used(edge))
				   {
					   tri.marked = 1;
					   my_hanging = true;
				   }
			   }
			 }
		 }
		 if (my_hanging)
		 {
			 hanging = true;
		 }
	});
		return hanging;
	  }



	  public static int MarkHangingQuads(Array<MarkedQuad> mquads, INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges)
	  {
		int hanging = 0;
		for (int i = 1; i <= mquads.Size(); i++)
		{
		if (mquads.Elem(i).marked)
		{
			hanging = 1;
			continue;
		}

		INDEX_2 edge1 = new INDEX_2(mquads.Get(i).pnums[0], mquads.Get(i).pnums[1]);
		INDEX_2 edge2 = new INDEX_2(mquads.Get(i).pnums[2], mquads.Get(i).pnums[3]);
		edge1.Sort();
		edge2.Sort();
		if (cutedges.Used(edge1) || cutedges.Used(edge2))
		{
			mquads.Elem(i).marked = 1;
				mquads.Elem(i).markededge = 0;
			hanging = 1;
				continue;
		}

			// he/sz: second case: split horizontally
			INDEX_2 edge3 = new INDEX_2(mquads.Get(i).pnums[1], mquads.Get(i).pnums[3]);
			INDEX_2 edge4 = new INDEX_2(mquads.Get(i).pnums[2], mquads.Get(i).pnums[0]);

			edge3.Sort();
			edge4.Sort();
			if (cutedges.Used(edge3) || cutedges.Used(edge4))
			{
			  mquads.Elem(i).marked = 1;
			  mquads.Elem(i).markededge = 1;
			  hanging = 1;
			  continue;
			}

		}
		return hanging;
	  }



	  public static void ConnectToNodeRec(int node, int tonode, TABLE<int> conto, Array<int> connecttonode)
	  {
		//  (*testout) << "connect " << node << " to " << tonode << endl;
		for (int i = 1; i <= conto.EntrySize(node); i++)
		{
		int n2 = conto.Get(node, i);
		if (!connecttonode.Get(n2))
		{
			connecttonode.Elem(n2) = tonode;
			ConnectToNodeRec(n2, tonode, conto, connecttonode);
		}
		}
	  }




	  public static Array<MarkedTet> mtets = new Array<MarkedTet>();
	  public static Array<MarkedPrism> mprisms = new Array<MarkedPrism>();
	  public static Array<MarkedIdentification> mids = new Array<MarkedIdentification>();
	  public static Array<MarkedTri> mtris = new Array<MarkedTri>();
	  public static Array<MarkedQuad> mquads = new Array<MarkedQuad>();


	  public static void WriteMarkedElements(ostream ost)
	  {
		ost << "Marked Elements\n";

		ost << mtets.Size() << "\n";
		for (int i = 0; i < mtets.Size(); i++)
		{
		  ost << mtets[i];
		}

		ost << mprisms.Size() << "\n";
		for (int i = 0; i < mprisms.Size(); i++)
		{
		  ost << mprisms[i];
		}

		ost << mids.Size() << "\n";
		for (int i = 0; i < mids.Size(); i++)
		{
		  ost << mids[i];
		}

		ost << mtris.Size() << "\n";
		for (int i = 0; i < mtris.Size(); i++)
		{
		  ost << mtris[i];
		}

		ost << mquads.Size() << "\n";
		for (int i = 0; i < mquads.Size(); i++)
		{
		  ost << mquads[i];
		}
		ost << "\n";
	  }

	  public static bool ReadMarkedElements(istream ist, Mesh mesh)
	  {
		string auxstring = "";
		if (ist != null)
		{
		  ist >> auxstring;
		}

		if (auxstring != "Marked")
		{
		  return false;
		}

		if (ist != null)
		{
		  ist >> auxstring;
		}

		if (auxstring != "Elements")
		{
		  return false;
		}

		int size;

		ist >> size;
		mtets.SetSize(size);
		for (int i = 0; i < size; i++)
		{
			ist >> mtets[i];
			if (mtets[i].pnums[0] > mesh.GetNV() || mtets[i].pnums[1] > mesh.GetNV() || mtets[i].pnums[2] > mesh.GetNV() || mtets[i].pnums[3] > mesh.GetNV())
			{
			  return false;
			}
		}

		ist >> size;
		mprisms.SetSize(size);
		for (int i = 0; i < size; i++)
		{
		  ist >> mprisms[i];
		}

		ist >> size;
		mids.SetSize(size);
		for (int i = 0; i < size; i++)
		{
		  ist >> mids[i];
		}

		ist >> size;
		mtris.SetSize(size);
		for (int i = 0; i < size; i++)
		{
		  ist >> mtris[i];
		}

		ist >> size;
		mquads.SetSize(size);
		for (int i = 0; i < size; i++)
		{
		  ist >> mquads[i];
		}

		return true;
	  }





	  public static void BisectTetsCopyMesh(Mesh mesh, NetgenGeometry UnnamedParameter, BisectionOptions opt, Array< Array<int,PointIndex.BASE> > idmaps, string refinfofile)
	  {
		// mtets.SetName ("bisection, tets");
		// mprisms.SetName ("bisection, prisms");
		// mtris.SetName ("bisection, trigs");
		// nmquads.SetName ("bisection, quads");
		// mids.SetName ("bisection, identifications");

		//int np = mesh.GetNP();
		int ne = mesh.GetNE();
		int nse = mesh.GetNSE();

		/*
		  if (mtets.Size() + mprisms.Size() == mesh.GetNE())
		  return;
		*/

		bool readok = false;

		if (refinfofile != "")
		{
		PrintMessage(3, "Reading marked-element information from \"", refinfofile, "\"");
		ifstream ist = new ifstream(refinfofile);

		readok = ReadMarkedElements(ist, mesh);

		ist.close();
		}

		if (!readok)
		{
		PrintMessage(3, "resetting marked-element information");
		mtets.SetSize(0);
		mprisms.SetSize(0);
		mids.SetSize(0);
		mtris.SetSize(0);
		mquads.SetSize(0);


		INDEX_2_HASHTABLE<int> shortedges = new INDEX_2_HASHTABLE<int>(100);
		for (int i = 1; i <= ne; i++)
		{
			Element el = mesh.VolumeElement(i);
			if (el.GetType() == ELEMENT_TYPE.PRISM || el.GetType() == ELEMENT_TYPE.PRISM12)
			{
			for (int j = 1; j <= 3; j++)
			{
				INDEX_2 se = new INDEX_2(el.PNum(j), el.PNum(j + 3));
				se.Sort();
				shortedges.Set(se, 1);
			}
			}
		}



		// INDEX_2_HASHTABLE<int> edgenumber(np);
		INDEX_2_CLOSED_HASHTABLE<int> edgenumber = new INDEX_2_CLOSED_HASHTABLE<int>((uint)(9 * ne+4 * nse));

		BTSortEdges(mesh, idmaps, edgenumber);


		for (int i = 1; i <= ne; i++)
		{
			Element el = mesh.VolumeElement(i);

			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TET:
			  case ELEMENT_TYPE.TET10:
			  {
			  // if tet has short edge, it is handled as degenerated prism

			  int foundse = 0;
			  for (int j = 1; j <= 3; j++)
			  {
				for (int k = j + 1; k <= 4; k++)
				{
				INDEX_2 se = new INDEX_2(el.PNum(j), el.PNum(k));
				se.Sort();
				if (shortedges.Used(se))
				{
					//		      cout << "tet converted to prism" << endl;

					foundse = 1;
					int p3 = 1;
					while (p3 == j || p3 == k)
					{
					  p3++;
					}
					int p4 = 10 - j - k - p3;

					// even permutation ?
					int[] pi = new int[4];
					pi[0] = j;
					pi[1] = k;
					pi[2] = p3;
					pi[3] = p4;
					int cnt = 0;
					for (int l = 1; l <= 4; l++)
					{
					  for (int m = 0; m < 3; m++)
					  {
					if (pi[m] > pi[m + 1])
					{
						Swap(ref pi[m], ref pi[m + 1]);
						cnt++;
					}
					  }
					}
					if ((cnt % 2) != 0)
					{
					  Swap(ref p3, ref p4);
					}

					Element hel = new Element(el);
					hel.PNum(1) = el.PNum(j);
					hel.PNum(2) = el.PNum(k);
					hel.PNum(3) = el.PNum(p3);
					hel.PNum(4) = el.PNum(p4);

					MarkedPrism mp = new MarkedPrism();
					BTDefineMarkedPrism(hel, edgenumber, mp);
					mp.matindex = el.GetIndex();
					mprisms.Append(mp);
				}
				}
			  }
			  if (foundse == 0)
			  {
				  MarkedTet mt = new MarkedTet();
				  BTDefineMarkedTet(el, edgenumber, mt);
				  mt.matindex = el.GetIndex();
				  mtets.Append(mt);
			  }
			  break;
			  }
			  case ELEMENT_TYPE.PYRAMID:
			  {
			  // eventually rotate
			  MarkedPrism mp = new MarkedPrism();

			  INDEX_2 se = new INDEX_2(el.PNum(1), el.PNum(2));
			  se.Sort();
			  if (shortedges.Used(se))
			  {
				  Element hel = new Element(el);
				  hel.PNum(1) = el.PNum(2);
				  hel.PNum(2) = el.PNum(3);
				  hel.PNum(3) = el.PNum(4);
				  hel.PNum(4) = el.PNum(1);
				  BTDefineMarkedPrism(hel, edgenumber, mp);
			  }
			  else
			  {
				  BTDefineMarkedPrism(el, edgenumber, mp);
			  }

			  mp.matindex = el.GetIndex();
			  mprisms.Append(mp);
			  break;
			  }
			  case ELEMENT_TYPE.PRISM:
			  case ELEMENT_TYPE.PRISM12:
			  {
			  MarkedPrism mp = new MarkedPrism();
			  BTDefineMarkedPrism(el, edgenumber, mp);
			  mp.matindex = el.GetIndex();
			  mprisms.Append(mp);
			  break;
			  }
				  default:
					throw new Exception("Bisect, element type not handled in switch, 4");
			}
		}

		for (int i = 1; i <= nse; i++)
		{
			Element2d el = mesh.SurfaceElement(i);
			if (el.GetType() == ELEMENT_TYPE.TRIG || el.GetType() == ELEMENT_TYPE.TRIG6)
			{
			MarkedTri mt = new MarkedTri();
			BTDefineMarkedTri(el, edgenumber, mt);
			mtris.Append(mt);
			}
			else
			{
			MarkedQuad mq = new MarkedQuad();
			BTDefineMarkedQuad(el, edgenumber, mq);
			mquads.Append(mq);
			}

				if (mesh.GetDimension() == 3)
				{
					MarkedIdentification mi = new MarkedIdentification();
					for (int j = 0; j < idmaps.Size(); j++)
					{
					  if (BTDefineMarkedId(el, edgenumber, idmaps[j], mi))
					  {
						mids.Append(mi);
					  }
					}
				}
		}
			if (mesh.GetDimension() == 2)
			{
				for (SegmentIndex j = 1; j <= mesh.GetNSeg(); j++)
				{
					var seg = mesh[j];
					foreach (var map in idmaps)
					{
						if (seg[0] > 0 && seg[1] > 0 && (*map)[seg[0]] && (*map)[seg[1]])
						{
							MarkedIdentification mi = new MarkedIdentification();
							mi.np = 2;
							mi.pnums[0] = seg[0];
							mi.pnums[1] = seg[1];
							mi.pnums[2] = (*map)[seg[0]];
							mi.pnums[3] = (*map)[seg[1]];
							var min1 = mi.pnums[0] < mi.pnums[1] ? mi.pnums[0] : mi.pnums[1];
							var min2 = mi.pnums[2] < mi.pnums[3] ? mi.pnums[2] : mi.pnums[3];
							if (min1 > min2)
							{
							  continue;
							}
							mi.marked = 0;
							mi.markededge = 0;
							mi.incorder = false;
							mids.Append(mi);
						}
					}
				}
			}
		}




		mesh.mlparentelement.SetSize(ne);
		for (int i = 1; i <= ne; i++)
		{
		  mesh.mlparentelement.Elem(i) = 0;
		}
		mesh.mlparentsurfaceelement.SetSize(nse);
		for (int i = 1; i <= nse; i++)
		{
		  mesh.mlparentsurfaceelement.Elem(i) = 0;
		}

		if (printmessage_importance > 0)
		{
		  ostringstream str1 = new ostringstream();
		  ostringstream str2 = new ostringstream();
		  str1 << "copied " << mtets.Size() << " tets, " << mprisms.Size() << " prisms";
		  str2 << "       " << mtris.Size() << " trigs, " << mquads.Size() << " quads";

		  PrintMessage(4, str1.str());
		  PrintMessage(4, str2.str());
		}
	  }


	  /*
	  void UpdateEdgeMarks2(Mesh & mesh,
				const Array< Array<int,PointIndex::BASE>* > & idmaps)
	  {
	    Array< Array<MarkedTet>*,PointIndex::BASE > mtets_old(mesh.GetNP());
	    Array< Array<MarkedPrism>*,PointIndex::BASE > mprisms_old(mesh.GetNP());
	    Array< Array<MarkedIdentification>*,PointIndex::BASE > mids_old(mesh.GetNP());
	    Array< Array<MarkedTri>*,PointIndex::BASE > mtris_old(mesh.GetNP());
	    Array< Array<MarkedQuad>*,PointIndex::BASE > mquads_old(mesh.GetNP());
	
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      mtets_old[i] = new Array<MarkedTet>;
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      mprisms_old[i] = new Array<MarkedPrism>;
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      mids_old[i] = new Array<MarkedIdentification>;
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      mtris_old[i] = new Array<MarkedTri>;
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      mquads_old[i] = new Array<MarkedQuad>;
	
	    for(int i=0; i<mtets.Size(); i++)
	      mtets_old[mtets[i].pnums[0]]->Append(mtets[i]);
	    for(int i=0; i<mprisms.Size(); i++)
	      mprisms_old[mprisms[i].pnums[0]]->Append(mprisms[i]);
	    for(int i=0; i<mids.Size(); i++)
	      mids_old[mids[i].pnums[0]]->Append(mids[i]);
	    for(int i=0; i<mtris.Size(); i++)
	      {
		(*testout) << "i " << i << endl;
		(*testout) << "mtris[i] " << mtris[i].pnums[0] << " " << mtris[i].pnums[1] << " " << mtris[i].pnums[2] << endl; 
		mtris_old[mtris[i].pnums[0]]->Append(mtris[i]);
	      }
	    for(int i=0; i<mquads.Size(); i++)
	      mquads_old[mquads[i].pnums[0]]->Append(mquads[i]);
	
	   
	    
	    int np = mesh.GetNP();
	    int ne = mesh.GetNE();
	    int nse = mesh.GetNSE();
	    int i, j, k, l, m;
	
	
	//       if (mtets.Size() + mprisms.Size() == mesh.GetNE())
	//       return;
	
	    
	
	    mtets.SetSize(0);
	    mprisms.SetSize(0);
	    mids.SetSize(0);
	    mtris.SetSize(0);
	    mquads.SetSize(0);
	
	
	    INDEX_2_HASHTABLE<int> shortedges(100);
	    for (i = 1; i <= ne; i++)
	      {
		const Element & el = mesh.VolumeElement(i);
		if (el.GetType() == PRISM ||
			el.GetType() == PRISM12)
		  {
			for (j = 1; j <= 3; j++)
			  {
			INDEX_2 se(el.PNum(j), el.PNum(j+3));
			se.Sort();
			shortedges.Set (se, 1);
			  }
		  }
	      }
	
	
	
	    // INDEX_2_HASHTABLE<int> edgenumber(np);
	    INDEX_2_CLOSED_HASHTABLE<int> edgenumber(9*ne+4*nse);  
	
	    BTSortEdges (mesh, idmaps, edgenumber);
	
	
	    for (i = 1; i <= ne; i++)
	      {
		const Element & el = mesh.VolumeElement(i);
	
		switch (el.GetType())
		  {
		  case TET:
		  case TET10:
			{
			  // if tet has short edge, it is handled as degenerated prism
	
			  int foundse = 0;
			  for (j = 1; j <= 3; j++)
			for (k = j+1; k <= 4; k++)
			  {
				INDEX_2 se(el.PNum(j), el.PNum(k));
				se.Sort();
				if (shortedges.Used (se))
				  {
	//		      cout << "tet converted to prism" << endl;
	
				foundse = 1;
				int p3 = 1;
				while (p3 == j || p3 == k)
				  p3++;
				int p4 = 10 - j - k - p3;
	
				// even permutation ?
				int pi[4];
				pi[0] = j;
				pi[1] = k;
				pi[2] = p3;
				pi[3] = p4;
				int cnt = 0;
				for (l = 1; l <= 4; l++)
				  for (m = 0; m < 3; m++)
					if (pi[m] > pi[m+1])
					  {
					Swap (pi[m], pi[m+1]);
					cnt++;
					  }
				if (cnt % 2)
				  Swap (p3, p4);
	
				Element hel = el;
				hel.PNum(1) = el.PNum(j);
				hel.PNum(2) = el.PNum(k);
				hel.PNum(3) = el.PNum(p3);
				hel.PNum(4) = el.PNum(p4);
	
				MarkedPrism mp;
	
				BTDefineMarkedPrism (hel, edgenumber, mp);
				mp.matindex = el.GetIndex();
				mprisms.Append (mp);
				  }
			  }
			  if (!foundse)
			{
			  MarkedTet mt;
	
			  int oldind = -1;
			  for(l = 0; oldind < 0 && l<mtets_old[el[0]]->Size(); l++)
				if(el[1] == (*mtets_old[el[0]])[l].pnums[1] &&
				   el[2] == (*mtets_old[el[0]])[l].pnums[2] &&
				   el[3] == (*mtets_old[el[0]])[l].pnums[3])
				  oldind = l;
	
			  if(oldind >= 0)
				mtets.Append((*mtets_old[el[0]])[oldind]);
			  else
				{
				  BTDefineMarkedTet (el, edgenumber, mt);
				  mt.matindex = el.GetIndex();
				  mtets.Append (mt);
				}
			}
			  break;
			}
		  case PYRAMID:
			{
			  // eventually rotate
			  MarkedPrism mp;
	
			  INDEX_2 se(el.PNum(1), el.PNum(2));
			  se.Sort();
			  if (shortedges.Used (se))
			{
			  Element hel = el;
			  hel.PNum(1) = el.PNum(2);
			  hel.PNum(2) = el.PNum(3);
			  hel.PNum(3) = el.PNum(4);
			  hel.PNum(4) = el.PNum(1);
			  BTDefineMarkedPrism (hel, edgenumber, mp);
			}
			  else
			{
			  BTDefineMarkedPrism (el, edgenumber, mp);
			}
	
			  mp.matindex = el.GetIndex();
			  mprisms.Append (mp);
			  break;
			}
		  case PRISM:
		  case PRISM12:
			{
			  MarkedPrism mp;
			  BTDefineMarkedPrism (el, edgenumber, mp);
			  mp.matindex = el.GetIndex();
			  mprisms.Append (mp);
			  break;
			}
		  }
	      }
	
	    for (i = 1; i <= nse; i++)
	      {
		const Element2d & el = mesh.SurfaceElement(i);
		if (el.GetType() == TRIG ||
			el.GetType() == TRIG6)
		  {
			MarkedTri mt;
			BTDefineMarkedTri (el, edgenumber, mt);
			mtris.Append (mt);
		  }
		else
		  {
			MarkedQuad mq;
			BTDefineMarkedQuad (el, edgenumber, mq);
			mquads.Append (mq);
		  }
	
		MarkedIdentification mi;
	
	
	
		for(j=0; j<idmaps.Size(); j++)
		  if(BTDefineMarkedId(el, edgenumber, *idmaps[j], mi))
			{
			  mids.Append(mi);
	
			  int oldind = -1;
			  for(l = 0; oldind < 0 && l<mids_old[mi.pnums[0]]->Size(); l++)
			{
			  bool equal = true;
			  for(int m = 1; equal && m < mi.np; m++)
				equal = (mi.pnums[m] == (*mids_old[el[0]])[l].pnums[m]);
			  if(equal)
				oldind = l;
			}
	
			  if(oldind >= 0)
			mids.Last() = (*mids_old[mi.pnums[0]])[oldind];
			}
	
	      }
	
	
	
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      delete mtets_old[i];
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      delete mprisms_old[i];
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      delete mids_old[i];
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      delete mtris_old[i];
	    for(int i=PointIndex::BASE; i<mesh.GetNP()+PointIndex::BASE; i++)
	      delete mquads_old[i];
	  }
	*/


	  public static void UpdateEdgeMarks(Mesh mesh, Array< Array<int,PointIndex.BASE> > idmaps)
	  //const Array < Array<Element>* > & elements_before,
	  //const Array < Array<int>* > & markedelts_num,
	  //		const Array < Array<Element2d>* > & surfelements_before,
	  //		const Array < Array<int>* > & markedsurfelts_num)
	  {
		/*
		T_MTETS mtets_old; mtets_old.Copy(mtets);
		T_MPRISMS mprisms_old; mprisms_old.Copy(mprisms);
		T_MIDS mids_old; mids_old.Copy(mids);
		T_MTRIS mtris_old; mtris_old.Copy(mtris);
		T_MQUADS mquads_old; mquads_old.Copy(mquads);
		*/
		Array<MarkedTet> mtets_old = new Array<MarkedTet>(mtets);
		Array<MarkedPrism> mprisms_old = new Array<MarkedPrism>(mprisms);
		Array<MarkedIdentification> mids_old = new Array<MarkedIdentification>(mids);
		Array<MarkedTri> mtris_old = new Array<MarkedTri>(mtris);
		Array<MarkedQuad> mquads_old = new Array<MarkedQuad>(mquads);



		mtets.SetSize(0);
		mprisms.SetSize(0);
		mids.SetSize(0);
		mtris.SetSize(0);
		mquads.SetSize(0);

		//int nv = mesh.GetNV();


		INDEX_2_CLOSED_HASHTABLE<int> edgenumber = new INDEX_2_CLOSED_HASHTABLE<int>(9 * mesh.GetNE() + 4 * mesh.GetNSE());

		int maxnum = BTSortEdges(mesh, idmaps, edgenumber);

		for (int m = 0; m < mtets_old.Size(); m++)
		{
		MarkedTet mt = mtets_old[m];

		//(*testout) << "old mt " << mt;

		INDEX_2 edge = new INDEX_2(mt.pnums[mt.tetedge1], mt.pnums[mt.tetedge2]);
		edge.Sort();
		if (edgenumber.Used(edge))
		{
			int val = edgenumber.Get(edge);
			//(*testout) << "set voledge " << edge << " from " << val;
			if (val <= maxnum)
			{
			val += 2 * maxnum;
			edgenumber.Set(edge, val);
			}
			else if (val <= 2 * maxnum)
			{
			val += maxnum;
			edgenumber.Set(edge, val);
			}
			//(*testout) << " to " << val << endl;
		}

		for (int k = 0; k < 4; k++)
		{
		  for (int i = 0; i < 3; i++)
		  {
			for (int j = i + 1; i != k && j < 4; j++)
			{
			  if (j != k && (int)mt.faceedges[k] == 6 - k - i - j)
			  {
			  edge[0] = mt.pnums[i];
			  edge[1] = mt.pnums[j];
			  edge.Sort();
			  if (edgenumber.Used(edge))
			  {
				  int val = edgenumber.Get(edge);
				  //(*testout) << "set faceedge " << edge << " from " << val;
				  if (val <= maxnum)
				  {
				  val += maxnum;
				  edgenumber.Set(edge, val);
				  }
				  //(*testout) << " to " << val << endl;
			  }
			  }
			}
		  }
		}
		}




		for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
		{
		Element el = mesh[ei];

		//int pos = elements_before[el[0]]->Pos(el);
		//int elnum = (pos >= 0) ? (*markedelts_num[el[0]])[pos] : -1;

		switch (el.GetType())
		{
		  case ELEMENT_TYPE.TET:
		  case ELEMENT_TYPE.TET10:
		  {
			  //if(elnum >= 0)
			  // {
			  //   mtets.Append(mtets_old[elnum]);
			  // } 
			  //else
			  // {
			  MarkedTet mt = new MarkedTet();
			  BTDefineMarkedTet(el, edgenumber, mt);
			  mt.matindex = el.GetIndex();

			  mtets.Append(mt);

			  //(*testout) << "mtet " << mtets.Last() << endl;
			  break;
		  }
		  case ELEMENT_TYPE.PYRAMID:
		  {
			  cerr << "Refinement :: UpdateEdgeMarks not yet implemented for pyramids" << "\n";
			  break;
		  }

		  case ELEMENT_TYPE.PRISM:
		  case ELEMENT_TYPE.PRISM12:
		  {
			  cerr << "Refinement :: UpdateEdgeMarks not yet implemented for prisms" << "\n";
			  break;
		  }
			  default:
				throw new Exception("Bisect, element type not handled in switch, 6");
		}

		}



		 for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
		 {
		 Element2d el = mesh[sei];

		 /*
		 for(int k=0; k<3; k++)
		   auxind3[k] = el[k];
	
		 auxind3.Sort();
		 
		 int pos = oldfaces[auxind3[0]]->Pos(auxind3);
		 if(pos < 0)
		   cout << "UIUIUI" << endl;
		 */	 

		 switch (el.GetType())
		 {
		   case ELEMENT_TYPE.TRIG:
		   case ELEMENT_TYPE.TRIG6:
		   {
			   MarkedTri mt = new MarkedTri();
			   BTDefineMarkedTri(el, edgenumber, mt);
			   mtris.Append(mt);
			   break;
		   }

		   case ELEMENT_TYPE.QUAD:
		   case ELEMENT_TYPE.QUAD6:
		   {
			   MarkedQuad mt = new MarkedQuad();
			   BTDefineMarkedQuad(el, edgenumber, mt);
			   mquads.Append(mt);
			   break;
		   }
			   default:
				 throw new Exception("Bisect, element type not handled in switch, 5");
		 }


		MarkedIdentification mi = new MarkedIdentification();
		for (int j = 0; j < idmaps.Size(); j++)
		{
		  if (BTDefineMarkedId(el, edgenumber, idmaps[j], mi))
		  {
			mids.Append(mi);
		  }
		}


		 /*
		 int pos = surfelements_before[el[0]]->Pos(el);
		 int elnum = (pos >= 0) ? (*markedsurfelts_num[el[0]])[pos] : -1;
		 
		 
		 switch (el.GetType())
		   {
		   case TRIG:
		   case TRIG6:
		     {
		       if(elnum >= 0)
			 mtris.Append(mtris_old[elnum]);
		       else
			 {
			   MarkedTri mt;
			   BTDefineMarkedTri (el, edgenumber, mt);
			   mtris.Append (mt);
			   (*testout) << "(new) ";
			 }
		       (*testout) << "mtri " << mtris.Last();
		       break;
		     }
		     
		   case QUAD:
		   case QUAD6:
		     {
		       if(elnum >= 0)
			 mquads.Append(mquads_old[elnum]);
		       else
			 {
			   MarkedQuad mt;
			   BTDefineMarkedQuad (el, edgenumber, mt);
			   mquads.Append (mt);
			 }
		       break;
		     }
		   }
		 */
		 }

		 /*
		 for(int i=0; i<oldfaces.Size(); i++)
		   {
		 delete oldfaces[i];
		 delete oldmarkededges[i];
		   }
		 */

	  }




	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer = NgProfiler.CreateTimer("Bisect");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer1 = NgProfiler.CreateTimer("Bisect 1");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer1a = NgProfiler.CreateTimer("Bisect 1a");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer1b = NgProfiler.CreateTimer("Bisect 1b");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer2 = NgProfiler.CreateTimer("Bisect 2");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer2a = NgProfiler.CreateTimer("Bisect 2a");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer2b = NgProfiler.CreateTimer("Bisect 2b");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer3 = NgProfiler.CreateTimer("Bisect 3");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer3a = NgProfiler.CreateTimer("Bisect 3a");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer3b = NgProfiler.CreateTimer("Bisect 3b");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer_bisecttet = NgProfiler.CreateTimer("Bisect tets");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer_bisecttrig = NgProfiler.CreateTimer("Bisect trigs");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_timer_bisectsegms = NgProfiler.CreateTimer("Bisect segms");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_localizetimer = NgProfiler.CreateTimer("localize edgepoints");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static bool Bisect_repaired_once;
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Bisect_reptimer = NgProfiler.CreateTimer("check/repair");





	/*
	   Philippose Rajan - 11 June 2009
	
	   Function to calculate the surface normal at a given 
	   vertex of a surface element, with respect to that 
	   surface element.
	
	   This function is used by the boundary layer generation 
	   function, in order to calculate the effective direction 
	   in which the prismatic layer should grow
	*/
	   public static void GetSurfaceNormal(Mesh mesh, Element2d el, int Vertex, ref Vec3d SurfaceNormal)
	   {
		  int Vertex_A;
		  int Vertex_B;

		  Vertex_A = Vertex + 1;
		  if (Vertex_A > el.GetNP())
		  {
			  Vertex_A = 1;
		  }

		  Vertex_B = Vertex - 1;
		  if (Vertex_B <= 0)
		  {
			  Vertex_B = el.GetNP();
		  }

		  Vec3d Vect_A = new Vec3d();
		  Vec3d Vect_B = new Vec3d();

		  Vect_A = mesh[el.PNum(Vertex_A)] - mesh[el.PNum(Vertex)];
		  Vect_B = mesh[el.PNum(Vertex_B)] - mesh[el.PNum(Vertex)];

		  SurfaceNormal = Cross(new netgen.Vec3d(Vect_A), new netgen.Vec3d(Vect_B));
		  SurfaceNormal.Normalize();
	   }

	  //   bool rational = true;


	  internal static void ComputeGaussRule(int n, Array<double> xi, Array<double> wi)
	  {
		xi.SetSize(n);
		wi.SetSize(n);

		int m = (n + 1) / 2;
		double p1;
		double p2;
		double p3;
		double pp;
		double z;
		double z1;
		for (int i = 1; i <= m; i++)
		{
		z = ngsimd.GlobalMembers.cos(DefineConstants.M_PI * (i - 0.25) / (n + 0.5));
		while (true)
		{
			p1 = 1;
			p2 = 0;
			for (int j = 1; j <= n; j++)
			{
			p3 = p2;
			p2 = p1;
			p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
			}
			// p1 is legendre polynomial

			pp = n * (z * p1 - p2) / (z * z - 1);
			z1 = z;
			z = z1 - p1 / pp;

			if (ngsimd.GlobalMembers.fabs(z - z1) < 1e-14)
			{
				break;
			}
		}

		xi[i - 1] = 0.5 * (1 - z);
		xi[n - i] = 0.5 * (1 + z);
		wi[i - 1] = wi[n - i] = 1.0 / ((1 - z * z) * pp * pp);
		}
	  }



	  // compute edge bubbles up to order n, x \in (-1, 1)
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  internal static void CalcEdgeShape<T>(int n, T x, T[] shape)
	  {
		T p1 = x;
		T p2 = -1;
		T p3 = 0;
		for (int j = 2; j <= n; j++)
		{
		p3 = p2;
		p2 = p1;
		p1 = ((2 * j - 3) * x * p2 - (j - 3) * p3) / j;
		shape[j - 2] = p1;
		}
	  }
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T, typename FUNC>
	  internal static void CalcEdgeShapeLambda<T, FUNC>(int n, T x, FUNC func)
	  {
		T p1 = new T(x);
		T p2 = new T(-1.0);
		T p3 = new T(0.0);
		for (int j = 2; j <= n; j++)
		{
		p3 = p2;
		p2 = p1;
		p1 = ((2 * j - 3) * x * p2 - (j - 3) * p3) / j;
		func(j - 2, p1);
		}
	  }


//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  internal static void CalcEdgeDx<T>(int n, T x, T[] dshape)
	  {
		T p1 = x;
		T p2 = -1;
		T p3 = 0;
		T p1dx = 1;
		T p2dx = 0;
		T p3dx = 0;

		for (int j = 2; j <= n; j++)
		{
		p3 = p2;
		p2 = p1;
		p3dx = p2dx;
		p2dx = p1dx;

		p1 = ((2 * j - 3) * x * p2 - (j - 3) * p3) / j;
		p1dx = ((2 * j - 3) * (x * p2dx + p2) - (j - 3) * p3dx) / j;

		dshape[j - 2] = p1dx;
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  internal static void CalcEdgeShapeDx<T>(int n, T x, T[] shape, T[] dshape)
	  {
		T p1 = x;
		T p2 = -1;
		T p3 = 0;
		T p1dx = 1;
		T p2dx = 0;
		T p3dx = 0;

		for (int j = 2; j <= n; j++)
		{
		p3 = p2;
		p2 = p1;
		p3dx = p2dx;
		p2dx = p1dx;

		p1 = ((2 * j - 3) * x * p2 - (j - 3) * p3) / j;
		p1dx = ((2 * j - 3) * (x * p2dx + p2) - (j - 3) * p3dx) / j;

		shape[j - 2] = p1;
		dshape[j - 2] = p1dx;
		}
	  }

	  // compute L_i(x/t) * t^i
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: There is no equivalent in C# to templates on variables:
	  internal static bool CalcScaledEdgeShape_init = false;
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static double[][] CalcScaledEdgeShape_coefs = RectangularArrays.RectangularDoubleArray(100, 2);

	  internal static void CalcScaledEdgeShape(int n, T x, T t, T[] shape)
	  {
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static bool init = false;
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static double coefs[100][2];
		if (!CalcScaledEdgeShape_init)
		{
			for (int j = 0; j < 100; j++)
			{
				CalcScaledEdgeShape_coefs[j][0] = (double)(2 * j + 1) / (j + 2);
				CalcScaledEdgeShape_coefs[j][1] = -(double)(j - 1) / (j + 2);
			}
			CalcScaledEdgeShape_init = true;
		}

		T p1 = x;
		T p2 = -1;
		T p3 = 0;
		T tt = t * t;
		for (int j = 0; j <= n - 2; j++)
		{
		p3 = p2;
		p2 = p1;
			p1 = CalcScaledEdgeShape_coefs[j][0] * x * p2 + CalcScaledEdgeShape_coefs[j][1] * tt * p3;
		// p1=( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
		shape[j] = p1;
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T, typename FUNC>
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: There is no equivalent in C# to templates on variables:
	  internal static bool CalcScaledEdgeShapeLambda_init = false;
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static double[][] CalcScaledEdgeShapeLambda_coefs = RectangularArrays.RectangularDoubleArray(100, 2);

	  internal static void CalcScaledEdgeShapeLambda(int n, T x, T t, FUNC func)
	  {
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static bool init = false;
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static double coefs[100][2];
		if (!CalcScaledEdgeShapeLambda_init)
		{
			for (int j = 0; j < 100; j++)
			{
				CalcScaledEdgeShapeLambda_coefs[j][0] = (double)(2 * j + 1) / (j + 2);
				CalcScaledEdgeShapeLambda_coefs[j][1] = -(double)(j - 1) / (j + 2);
			}
			CalcScaledEdgeShapeLambda_init = true;
		}

		T p1 = new T(x);
		T p2 = new T(-1.0);
		T p3 = new T(0.0);
		T tt = t * t;
		for (int j = 0; j <= n - 2; j++)
		{
		p3 = p2;
		p2 = p1;
			p1 = CalcScaledEdgeShapeLambda_coefs[j][0] * x * p2 + CalcScaledEdgeShapeLambda_coefs[j][1] * tt * p3;
		// p1=( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
		func(j, p1);
		}
	  }



//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <int DIST, typename T>
	  internal static void CalcScaledEdgeShapeDxDt<int DIST, T>(int n, T x, T t, T[] dshape)
	  {
		T p1 = x;
		T p2 = -1;
		T p3 = 0;
		T p1dx = 1;
		T p1dt = 0;
		T p2dx = 0;
		T p2dt = 0;
		T p3dx = 0;
		T p3dt = 0;

		for (int j = 0; j <= n - 2; j++)
		{
		p3 = p2;
		p3dx = p2dx;
		p3dt = p2dt;
		p2 = p1;
		p2dx = p1dx;
		p2dt = p1dt;

		p1 = ((2 * j + 1) * x * p2 - t * t * (j - 1) * p3) / (j + 2);
		p1dx = ((2 * j + 1) * (x * p2dx + p2) - t * t * (j - 1) * p3dx) / (j + 2);
		p1dt = ((2 * j + 1) * x * p2dt - (j - 1) * (t * t * p3dt + 2 * t * p3)) / (j + 2);

		// shape[j] = p1;
		dshape[DIST * j] = p1dx;
		dshape[DIST * j + 1] = p1dt;
		}
	  }


//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Tx, class Tres>
	  internal static void LegendrePolynomial<Tx, Tres>(int n, Tx x, Tres[] values)
	  {
		switch (n)
		{
		  case 0:
		values[0] = 1;
		break;
		  case 1:
		values[0] = 1;
		values[1] = x;
		break;

		  default:

		if (n < 0)
		{
			return;
		}

		Tx p1 = 1.0;
		Tx p2 = 0.0;
		Tx p3 = new default(Tx);

		values[0] = 1.0;
		for (int j = 1; j <= n; j++)
		{
			p3 = p2;
			p2 = p1;
			p1 = ((2.0 * j - 1.0) * x * p2 - (j - 1.0) * p3) / j;
			values[j] = p1;
		}
	  break;
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Tx, class Tt, class Tres>
	  internal static void ScaledLegendrePolynomial<Tx, Tt, Tres>(int n, Tx x, Tt t, Tres[] values)
	  {
		switch (n)
		{
		  case 0:
		values[0] = 1.0;
		break;

		  case 1:
		values[0] = 1.0;
		values[1] = x;
		break;

		  default:

		if (n < 0)
		{
			return;
		}

		Tx p1 = 1.0;
		Tx p2 = 0.0;
		Tx p3 = new default(Tx);
		values[0] = 1.0;
		for (int j = 1; j <= n; j++)
		{
			p3 = p2;
			p2 = p1;
			p1 = ((2.0 * j - 1.0) * x * p2 - t * t * (j - 1.0) * p3) / j;
			values[j] = p1;
		}
	  break;
		}
	  }
	  public static void JacobiPolynomial<S, T>(int n, S x, double alpha, double beta, T[] values)
	  {
		S p1 = 1.0;
		S p2 = 0.0;
		S p3 = new default(S);

		if (n >= 0)
		{
		  p2 = values[0] = 1.0;
		}
		if (n >= 1)
		{
		  p1 = values[1] = 0.5 * (2 * (alpha + 1) + (alpha + beta + 2) * (x - 1));
		}

		for (int i = 1; i < n; i++)
		{
			p3 = p2;
			p2 = p1;
			p1 = 1.0 / (2 * (i + 1) * (i + alpha + beta + 1) * (2 * i + alpha + beta)) * (((2 * i + alpha + beta + 1) * (alpha * alpha - beta * beta) + (2 * i + alpha + beta) * (2 * i + alpha + beta + 1) * (2 * i + alpha + beta + 2) * x) * p2 - 2 * (i + alpha) * (i + beta) * (2 * i + alpha + beta + 2) * p3);
			values[i + 1] = p1;
		}
	  }




//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class S, class St, class T>
	  public static void ScaledJacobiPolynomial<S, St, T>(int n, S x, St t, double alpha, double beta, T[] values)
	  {
		/*
		  S p1 = 1.0, p2 = 0.0, p3;
	
		  if (n >= 0) values[0] = 1.0;
		*/

		S p1 = new S(1.0);
		S p2 = new S(0.0);
		S p3 = new default(S);

		if (n >= 0)
		{
		  p2 = values[0] = 1.0;
		}
		if (n >= 1)
		{
		  p1 = values[1] = 0.5 * (2 * (alpha + 1) * t + (alpha + beta + 2) * (x - t));
		}

		for (int i = 1; i < n; i++)
		{
			p3 = p2;
			p2 = p1;
			p1 = 1.0 / (2 * (i + 1) * (i + alpha + beta + 1) * (2 * i + alpha + beta)) * (((2 * i + alpha + beta + 1) * (alpha * alpha - beta * beta) * t + (2 * i + alpha + beta) * (2 * i + alpha + beta + 1) * (2 * i + alpha + beta + 2) * x) * p2 - 2 * (i + alpha) * (i + beta) * (2 * i + alpha + beta + 2) * t * t * p3);
			values[i + 1] = p1;
		}
	  }


	  internal static Array<RecPol> jacpols2 = new Array<RecPol>();


	  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Tx, class Ty, class Ts>
	  internal static void CalcTrigShape<Tx, Ty, Ts>(int n, Tx x, Ty y, Ts[] shape)
	  {
		// cout << "calc trig shape" << endl;
		if (n < 3)
		{
			return;
		}
		Tx[] hx = Arrays.InitializeWithDefaultInstances<Tx>(50);
		Tx[] hy = Arrays.InitializeWithDefaultInstances<Tx>(50 * 50);

		jacpols2[2].EvaluateScaled(n - 3, x, 1 - y, hx);

		for (int ix = 0; ix <= n - 3; ix++)
		{
		  jacpols2[2 * ix + 5].Evaluate(n - 3, 2 * y - 1, hy + 50 * ix);
		}

		int ii = 0;

		Tx bub = (1 + x - y) * y * (1 - x - y);
		for (int ix = 0; ix <= n - 3; ix++)
		{
		  hx[ix] *= bub;
		}

		/*
		for (int iy = 0; iy <= n-3; iy++)
		  for (int ix = 0; ix <= n-3-iy; ix++)
		shape[ii++] = hx[ix]*hy[iy+50*ix];
		*/
		// change loops:
		for (int ix = 0; ix <= n - 3; ix++)
		{
		  for (int iy = 0; iy <= n - 3 - ix; iy++)
		  {
		shape[ii++] = hx[ix] * hy[iy + 50 * ix];
		  }
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  internal static void CalcTrigShapeDxDy<T>(int n, T x, T y, T[] dshape)
	  {
		if (n < 3)
		{
			return;
		}

		AutoDiff<2,T> adx = new AutoDiff<2,T>(x, 0);
		AutoDiff<2,T> ady = new AutoDiff<2,T>(y, 1);
		AutoDiff<2,T>[] res = Arrays.InitializeWithDefaultInstances<AutoDiff>(2000);
		CalcTrigShape(n, new AutoDiff<2,T>(adx), new AutoDiff<2,T>(ady), res[0]);
		int ndof = (n - 1) * (n - 2) / 2;
		for (int i = 0; i < ndof; i++)
		{
		dshape[2 * i] = res[i].DValue(0);
		dshape[2 * i + 1] = res[i].DValue(1);
		}
	  }


	  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Tx, class Ty, class Tt, class Tr>
	  internal static void CalcScaledTrigShape<Tx, Ty, Tt, Tr>(int n, Tx x, Ty y, Tt t, Tr[] shape)
	  {
		if (n < 3)
		{
			return;
		}

		Tx[] hx = Arrays.InitializeWithDefaultInstances<Tx>(50);
		Tx[] hy = Arrays.InitializeWithDefaultInstances<Tx>(50);
		ScaledJacobiPolynomial(n - 3, x, t - y, 2, 2, hx);

		int ii = 0;
		Tx bub = (t + x - y) * y * (t - x - y);
		for (int ix = 0; ix <= n - 3; ix++)
		{
			jacpols2[2 * ix + 5].EvaluateScaled(n - 3, 2 * y - 1, t, hy);
			for (int iy = 0; iy <= n - 3 - ix; iy++)
			{
			  shape[ii++] = bub * hx[ix] * hy[iy];
			}
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Tx, class Ty, class Tt, typename FUNC>
	  internal static void CalcScaledTrigShapeLambda<Tx, Ty, Tt, FUNC>(int n, Tx x, Ty y, Tt t, FUNC func)
	  {
		if (n < 3)
		{
			return;
		}
		int ii = 0;
		Tx bub = (t + x - y) * y * (t - x - y);
	jacpols2[2].EvaluateScaledLambda(n - 3, x, t - y, (int ix, Tx valx) =>
	{
		 jacpols2[2 * ix + 5].EvaluateScaledLambda(n - 3 - ix, 2 * y - 1, t, (int iy, Ty valy) =>
		 {
													 func(ii++, bub * valx * valy);
		 });
	});
	  }


	  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  internal static void CalcScaledTrigShapeDxDyDt<T>(int n, T x, T y, T t, T[] dshape)
	  {
		/*
		if (n < 3) return;
		AutoDiff<3,T> adx(x, 0);
		AutoDiff<3,T> ady(y, 1);
		AutoDiff<3,T> adt(t, 2);
		AutoDiff<3,T> res[2000];
		CalcScaledTrigShape (n, adx, ady, adt, &res[0]);
		int ndof = (n-1)*(n-2)/2;
		for (int i = 0; i < ndof; i++)
		  {
		dshape[3*i] = res[i].DValue(0);
		dshape[3*i+1] = res[i].DValue(1);
		dshape[3*i+2] = res[i].DValue(2);
		  }
		*/
		if (n < 3)
		{
			return;
		}
		AutoDiff<3,T> adx = new AutoDiff<3,T>(x, 0);
		AutoDiff<3,T> ady = new AutoDiff<3,T>(y, 1);
		AutoDiff<3,T> adt = new AutoDiff<3,T>(t, 2);
	CalcScaledTrigShapeLambda(n, new AutoDiff<3,T>(adx), new AutoDiff<3,T>(ady), new AutoDiff<3,T>(adt), (int i, AutoDiff<3,T> shape) =>
	{
								 dshape[3 * i] = shape.DValue(0);
								 dshape[3 * i + 1] = shape.DValue(1);
								 dshape[3 * i + 2] = shape.DValue(2);
	});
	  }


//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::GetCoefficients<2> (SurfaceElementInfo & info, Array<Vec<2>> & coefs) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::GetCoefficients (SurfaceElementInfo info, Array<Vec<2>> coefs);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::GetCoefficients<3> (SurfaceElementInfo & info, Array<Vec<3>> & coefs) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::GetCoefficients (SurfaceElementInfo info, Array<Vec<3>> coefs);


//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSegmentTransformation<2> (SegmentIndex elnr, int npts, const double * xi, uint sxi, double * x, uint sx, double * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSegmentTransformation (SegmentIndex elnr, int npts, double xi, uint sxi, ref double x, uint sx, ref double dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSegmentTransformation<3> (SegmentIndex elnr, int npts, const double * xi, uint sxi, double * x, uint sx, double * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSegmentTransformation (SegmentIndex elnr, int npts, double xi, uint sxi, ref double x, uint sx, ref double dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSegmentTransformation<2> (SegmentIndex elnr, int npts, const SIMD<double> * xi, uint sxi, SIMD<double> * x, uint sx, SIMD<double> * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSegmentTransformation (SegmentIndex elnr, int npts, SIMD<double> xi, uint sxi, SIMD<double> x, uint sx, SIMD<double> dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSegmentTransformation<3> (SegmentIndex elnr, int npts, const SIMD<double> * xi, uint sxi, SIMD<double> * x, uint sx, SIMD<double> * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSegmentTransformation (SegmentIndex elnr, int npts, SIMD<double> xi, uint sxi, SIMD<double> x, uint sx, SIMD<double> dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcSegmentTransformation<double> (double xi, SegmentIndex elnr, Point<3,double> * x, Vec<3,double> * dxdxi, bool * curved);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcSegmentTransformation (double xi, SegmentIndex elnr, Point<3,double> x, Vec<3,double> dxdxi, ref bool curved);



//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSurfaceTransformation<2> (SurfaceElementIndex elnr, int npts, const double * xi, uint sxi, double * x, uint sx, double * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int npts, double xi, uint sxi, ref double x, uint sx, ref double dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSurfaceTransformation<3> (SurfaceElementIndex elnr, int npts, const double * xi, uint sxi, double * x, uint sx, double * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int npts, double xi, uint sxi, ref double x, uint sx, ref double dxdxi, uint sdxdxi);


//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSurfaceTransformation<2> (SurfaceElementIndex elnr, int npts, const SIMD<double> * xi, uint sxi, SIMD<double> * x, uint sx, SIMD<double> * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int npts, SIMD<double> xi, uint sxi, SIMD<double> x, uint sx, SIMD<double> dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void CurvedElements::CalcMultiPointSurfaceTransformation<3> (SurfaceElementIndex elnr, int npts, const SIMD<double> * xi, uint sxi, SIMD<double> * x, uint sx, SIMD<double> * dxdxi, uint sdxdxi);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int npts, SIMD<double> xi, uint sxi, SIMD<double> x, uint sx, SIMD<double> dxdxi, uint sdxdxi);


//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointElementTransformation(ElementIndex elnr, int n, double xi, uint sxi, ref double x, uint sx, ref double dxdxi, uint sdxdxi);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void CurvedElements::CalcMultiPointElementTransformation(ElementIndex elnr, int n, SIMD<double> xi, uint sxi, SIMD<double> x, uint sx, SIMD<double> dxdxi, uint sdxdxi);

	  internal int[][] deltetfaces =
	  {
		  new int[] {1, 2, 3},
		  new int[] {2, 0, 3},
		  new int[] {0, 1, 3},
		  new int[] {1, 0, 2}
	  };





	  public static void AddDelaunayPoint(PointIndex newpi, Point3d newp, Array<DelaunayTet> tempels, Mesh mesh, BoxTree < 3> tettree, MeshNB meshnb, Array<Point < 3>> centers, Array<double> radi2, Array<int> connected, Array<int> treesearch, Array<int> freelist, SphereList list, IndexSet insphere, IndexSet closesphere)
	  {
		// static Timer t("Meshing3::AddDelaunayPoint"); RegionTimer reg(t);
		/*
		  find any sphere, such that newp is contained in
		*/
		DelaunayTet el = new DelaunayTet();
		int cfelind = -1;

		const Point < 3> * pp[4];
		Point < 3> pc;
		Point3d tpmin = new Point3d();
		Point3d tpmax = new Point3d();

		tettree.GetIntersecting(newp, newp, treesearch);

		double quot;
		double minquot = 1e20;

		foreach (var jjj in treesearch)
		{
		quot = Dist2(centers.Get(jjj), newp) / radi2.Get(jjj);

		if ((cfelind == -1 || quot < 0.99 * minquot) && quot < 1)
		{
			minquot = quot;
			el = tempels.Get(jjj);
			cfelind = jjj;
			if (minquot < 0.917632)
			{
			  break;
			}
		}
		}


		if (cfelind == -1)
		{
		PrintWarning("Delaunay, point not in any sphere");
		return;
		}

		/*
		  insphere:     point is in sphere -> delete element
		  closesphere:  point is close to sphere -> considered for same center
		*/

		// save overestimate
		insphere.SetMaxIndex(2 * tempels.Size() + 5 * mesh.GetNP());
		closesphere.SetMaxIndex(2 * tempels.Size() + 5 * mesh.GetNP());

		insphere.Clear();
		closesphere.Clear();


		insphere.Add(cfelind);

		int changed = 1;
		int nstarti = 1;
		int starti;

		while (changed != 0)
		{
		changed = 0;
		starti = nstarti;
		nstarti = insphere.GetArray().Size() + 1;


		// if point in sphere, then it is also closesphere
		for (int j = starti; j < nstarti; j++)
		{
			int helind = insphere.GetArray().Get(j);
			if (closesphere.IsIn(helind) == 0)
			{
			  closesphere.Add(helind);
			}
		}

		// add connected spheres to insphere - list
		for (int j = starti; j < nstarti; j++)
		{
			list.GetList(insphere.GetArray().Get(j), connected);
			for (int k = 0; k < connected.Size(); k++)
			{
			int celind = connected[k];

			if (tempels.Get(celind)[0] != -1 && insphere.IsIn(celind) == 0)
			{
				changed = 1;
				insphere.Add(celind);
			}
			}
		}

		// check neighbour-tets
		for (int j = starti; j < nstarti; j++)
		{
		  for (int k = 0; k < 4; k++)
		  {
			  int helind = insphere.GetArray().Get(j);
			  int nbind = meshnb.GetNB(helind, k);

			  if (nbind != 0 && insphere.IsIn(nbind) == 0)
			  {
			  //changed
			  //int prec = testout->precision();
			  //testout->precision(12);
			  //(*testout) << "val1 " << Dist2 (centers.Get(nbind), newp)
			  //	     << " val2 " << radi2.Get(nbind) * (1+1e-8)
			  //	     << " val3 " << radi2.Get(nbind)
			  //	     << " val1 / val3 " << Dist2 (centers.Get(nbind), newp)/radi2.Get(nbind) << endl;
			  //testout->precision(prec);
			  if (Dist2(centers.Get(nbind), newp) < radi2.Get(nbind) * (1 + 1e-8))
			  {
				closesphere.Add(nbind);
			  }

			  if (Dist2(centers.Get(nbind), newp) < radi2.Get(nbind) * (1 + 1e-12))
			  {
				  // point is in sphere -> remove tet
				  insphere.Add(nbind);
				  changed = 1;
			  }
			  else
			  {
				  INDEX_3 i3 = tempels.Get(helind).GetFace(k);

				  const Point < 3> & p1 = new mesh.Point(new PointIndex(i3.I1()));
				  const Point < 3> & p2 = new mesh.Point(new PointIndex(i3.I2()));
				  const Point < 3> & p3 = new mesh.Point(new PointIndex(i3.I3()));

				  Vec < 3> v1 = p2 - p1;
				  Vec < 3> v2 = p3 - p1;
				  Vec < 3> n = Cross(v1, v2);
				  n /= n.Length();

				  if (n * new Vec3d(p1, new mesh.Point(tempels.Get helind[k])) > 0)
				  {
				n *= -1;
				  }

				  double dist = n * new Vec3d(p1, newp);

				  if (dist > -1e-10) // 1e-10
				  {
				  insphere.Add(nbind);
				  changed = 1;
				  }
			  }
			  }
		  }
		}
		} // while (changed)

		// Array<Element> newels;
		Array<DelaunayTet> newels = new Array<DelaunayTet>();

		Element2d face = new Element2d(ELEMENT_TYPE.TRIG);

		foreach (int celind in insphere.GetArray())
		{
		  for (int k = 0; k < 4; k++)
		  {
		  int nbind = meshnb.GetNB(celind, k);

		  if (nbind == 0 || insphere.IsIn(nbind) == 0)
		  {
			  tempels.Get(celind).GetFace(k, face);

			  // Element newel(TET);
				  DelaunayTet newel = new DelaunayTet();
			  for (int l = 0; l < 3; l++)
			  {
					newel[l] = face[l];
			  }
				  newel[3] = newpi;

			  newels.Append(newel);

				  Vec < 3> v1 = mesh[face[1]] - mesh[face[0]];
				  Vec < 3> v2 = mesh[face[2]] - mesh[face[0]];
			  Vec < 3> n = Cross(v1, v2);

			  n.Normalize();
			  if (n * new Vec3d(new mesh.Point(face[0]), new mesh.Point(tempels.Get celind[k])) > 0)
			  {
			n *= -1;
			  }

				  double hval = n * (newp - mesh[face[0]]);

			  if (hval > -1e-12)
			  {
			  cerr << "vec to outer" << "\n";
			  (*testout) << "vec to outer, hval = " << hval << "\n";
			  (*testout) << "v1 x v2 = " << Cross(v1, v2) << "\n";
			  (*testout) << "facep: " << new mesh.Point(face[0]) << " " << new mesh.Point(face[1]) << " " << new mesh.Point(face[2]) << "\n";
			  }
		  }
		  }
		}

		meshnb.ResetFaceHT(10 * insphere.GetArray().Size() + 1);

		foreach (var celind in insphere.GetArray())
		{
		meshnb.Delete(celind);
		list.DeleteElement(celind);

		for (int k = 0; k < 4; k++)
		{
		  tempels.Elem[] celind = -1;
		}

			tettree.DeleteElement(celind);
		freelist.Append(celind);
		}



		bool hasclose = false;
		foreach (int ind in closesphere.GetArray())
		{
		if (insphere.IsIn(ind) == 0 && ngsimd.GlobalMembers.fabs(Dist2(centers.Get(ind), newp) - radi2.Get(ind)) < 1e-8)
		{
		  hasclose = true;
		}
		}

		for (int j = 1; j <= newels.Size(); j++)
		{
			var newel = newels.Get(j);
		int nelind;

		if (!freelist.Size())
		{
			tempels.Append(newel);
			nelind = tempels.Size();
		}
		else
		{
			nelind = freelist.Last();
			freelist.DeleteLast();

			tempels.Elem(nelind) = newel;
		}

		meshnb.Add(nelind);
		list.AddElement(nelind);

		for (int k = 0; k < 4; k++)
		{
		  pp[k] = &new mesh.Point(newel[k]);
		}

		if (CalcSphereCenter(pp[0], pc))
		{
			PrintSysError("Delaunay: New tet is flat");

			(*testout) << "new tet is flat" << "\n";
			for (int k = 0; k < 4; k++)
			{
			  (*testout) << newel[k] << " ";
			}
			(*testout) << "\n";
			for (int k = 0; k < 4; k++)
			{
			  (*testout) << *pp[k - 1] << " ";
			}
			(*testout) << "\n";
		}

		double r2 = Dist2(pp[0], pc);
		if (hasclose)
		{
		  for (int k = 1; k <= closesphere.GetArray().Size(); k++)
		  {
			  int csameind = closesphere.GetArray().Get(k);
			  if (insphere.IsIn(csameind) == 0 && ngsimd.GlobalMembers.fabs(r2 - radi2.Get(csameind)) < 1e-10 && Dist(pc, centers.Get(csameind)) < 1e-10)
			  {
			  pc = centers.Get(csameind);
			  r2 = radi2.Get(csameind);
			  list.ConnectElement(nelind, csameind);
			  break;
			  }
		  }
		}

		if (centers.Size() < nelind)
		{
			centers.Append(pc);
			radi2.Append(r2);
		}
		else
		{
			centers.Elem(nelind) = pc;
			radi2.Elem(nelind) = r2;
		}

		closesphere.Add(nelind);

		tpmax = tpmin = *pp[0];
		for (int k = 1; k <= 3; k++)
		{
			tpmin.SetToMin(pp[k]);
			tpmax.SetToMax(pp[k]);
		}
		tpmax.CopyFrom(tpmax + 0.01 * (tpmax - tpmin));
		tettree.Insert(tpmin, tpmax, nelind);
		}
	  }






	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  static Timer Delaunay1_t("Meshing3::Delaunay1");

	  public static void Delaunay1(Mesh mesh, MeshingParameters mp, AdFront3 adfront, Array<DelaunayTet> tempels, int oldnp, DelaunayTet startel, ref Point3d pmin, ref Point3d pmax)
	  {
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static Timer t("Meshing3::Delaunay1");
		RegionTimer reg = new RegionTimer(Delaunay1_t);

		Array<Point < 3>> centers = new Array<Point < 3>>();
		Array<double> radi2 = new Array<double>();

		Box < 3> bbox(Box < 3>.EMPTY_BOX);

		foreach (var face in adfront.Faces())
		{
		  foreach (PointIndex pi in face.Face().PNums())
		  {
			bbox.Add(new mesh.Point(pi));
		  }
		}

		foreach (PointIndex pi in mesh.LockedPoints())
		{
		  bbox.Add(new mesh.Point(pi));
		}

		pmin = bbox.PMin();
		pmax = bbox.PMax();


		Vec < 3> vdiag = pmax - pmin;
		// double r1 = vdiag.Length();
		double r1 = ngsimd.GlobalMembers.sqrt(3.0) * max3(vdiag(0), vdiag(1), vdiag(2));
		vdiag = Vec < 3> (r1, r1, r1);
		//double r2;

		Point < 3> pmin2 = pmin - 8 * vdiag;
		Point < 3> pmax2 = pmax + 8 * vdiag;

		Point < 3> cp1(pmin2), cp2(pmax2), cp3(pmax2), cp4(pmax2);
		cp2(0) = pmin2(0);
		cp3(1) = pmin2(1);
		cp4(2) = pmin2(2);


		uint np = mesh.GetNP();

		startel[0] = mesh.AddPoint(cp1);
		startel[1] = mesh.AddPoint(cp2);
		startel[2] = mesh.AddPoint(cp3);
		startel[3] = mesh.AddPoint(cp4);

		// flag points to use for Delaunay:
		BitArrayChar<PointIndex.BASE> usep = new BitArrayChar<PointIndex.BASE>(np);
		usep.Clear();

		foreach (var face in adfront.Faces())
		{
		  foreach (PointIndex pi in face.Face().PNums())
		  {
			usep.Set(pi);
		  }
		}

		for (uint i = oldnp + PointIndex.BASE; i < np + PointIndex.BASE; i++)
		{
		  usep.Set(i);
		}

		foreach (PointIndex pi in mesh.LockedPoints())
		{
		  usep.Set(pi);
		}


		Array<int> freelist = new Array<int>();

		int cntp = 0;

		MeshNB meshnb = new MeshNB(tempels, mesh.GetNP() + 5);
		SphereList list = new SphereList();

		pmin2 = pmin2 + 0.1 * (pmin2 - pmax2);
		pmax2 = pmax2 + 0.1 * (pmax2 - pmin2);

		BoxTree < 3> tettree(pmin2, pmax2);


		tempels.Append(startel);
		meshnb.Add(1);
		list.AddElement(1);
		Array<int> connected = new Array<int>();
		Array<int> treesearch = new Array<int>();

		Box < 3> tbox(Box < 3>.EMPTY_BOX);
		for (uint k = 0; k < 4; k++)
		{
		  tbox.Add(new mesh.Point(startel[k]));
		}
		Point < 3> tpmin = tbox.PMin();
		Point < 3> tpmax = tbox.PMax();

		tpmax = tpmax + 0.01 * (tpmax - tpmin);
		tettree.Insert(tpmin, tpmax, 1);

		Point < 3> pc;

		const Point < 3> * pp[4];
		for (int k = 0; k < 4; k++)
		{
		  pp[k] = &new mesh.Point(startel[k]);
		}
		CalcSphereCenter(pp[0], pc);

		centers.Append(pc);
		radi2.Append(Dist2(pp[0], pc));


		IndexSet insphere = new IndexSet(mesh.GetNP());
		IndexSet closesphere = new IndexSet(mesh.GetNP());

		// "random" reordering of points  (speeds a factor 3 - 5 !!!)
		Array<PointIndex, PointIndex.BASE, PointIndex> mixed = new Array<PointIndex, PointIndex.BASE, PointIndex>(np);
		int[] prims = {11, 13, 17, 19, 23, 29, 31, 37};
		int prim;

		{
		  int i = 0;
		  while (np % prims[i] == 0)
		  {
			  i++;
		  }
		  prim = prims[i];
		}

		for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End() - 4; pi++)
		{
		  mixed[pi] = new PointIndex((prim * pi) % np + PointIndex.BASE);
		}

		for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End() - 4; pi++)
		{
		if (pi % 1000 == 0 != null)
		{
			if (pi % 10000 == 0 != null)
			{
			  PrintDot('+');
			}
			else
			{
			  PrintDot('.');
			}
		}

		multithread.percent = 100.0 * pi / np;
		if (multithread.terminate)
		{
		  break;
		}

		PointIndex newpi = mixed[pi];

		if (usep.Test(newpi) == 0)
		{
		  continue;
		}

		cntp++;

		MeshPoint newp = mesh[newpi];

		AddDelaunayPoint(new netgen.PointIndex(newpi), newp, tempels, mesh, tettree, meshnb, centers, radi2, connected, treesearch, freelist, list, insphere, closesphere);

		}

		for (int i = tempels.Size(); i >= 1; i--)
		{
		  if (tempels.Get(i)[0] <= 0)
		  {
		tempels.DeleteElement(i);
		  }
		}

		PrintDot('\n');

		PrintMessage(3, "Points: ", cntp);
		PrintMessage(3, "Elements: ", tempels.Size());
		//   (*mycout) << cntp << " / " << tempels.Size() << " points/elements" << endl;

		/*
		  cout << "tempels: ";
		  tempels.PrintMemInfo(cout);
		  cout << "Searchtree: ";
		  tettree.Tree().PrintMemInfo(cout);
		  cout << "MeshNB: ";
		  meshnb.PrintMemInfo(cout);
		*/
	  }

	  public static ostream operator << (ostream ost, DelaunayTrig trig)
	  {
		ost << trig[0] << "-" << trig[1] << "-" << trig[2] << "\n";
		return ost;
	  }


	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int BlockFillLocalH_timer = NgProfiler.CreateTimer("Meshing2::BlockFill");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int BlockFillLocalH_timer1 = NgProfiler.CreateTimer("Meshing2::BlockFill 1");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int BlockFillLocalH_timer2 = NgProfiler.CreateTimer("Meshing2::BlockFill 2");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int BlockFillLocalH_timer3 = NgProfiler.CreateTimer("Meshing2::BlockFill 3");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int BlockFillLocalH_timer4 = NgProfiler.CreateTimer("Meshing2::BlockFill 4");




	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timer = NgProfiler.CreateTimer("Meshing2::Delaunay - total");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timerstart = NgProfiler.CreateTimer("Meshing2::Delaunay - start");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timerfinish = NgProfiler.CreateTimer("Meshing2::Delaunay - finish");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timer1 = NgProfiler.CreateTimer("Meshing2::Delaunay - incremental");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timer1a = NgProfiler.CreateTimer("Meshing2::Delaunay - incremental a");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timer1b = NgProfiler.CreateTimer("Meshing2::Delaunay - incremental b");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timer1c = NgProfiler.CreateTimer("Meshing2::Delaunay - incremental c");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int Delaunay_timer1d = NgProfiler.CreateTimer("Meshing2::Delaunay - incremental d");

	  ///
	  public static double GetTime()
	  {
		return (double)(clock() - starttimea) / DefineConstants.CLOCKS_PER_SEC;
	  }

	  public static void ResetTime()
	  {
		starttimea = clock();
	  }

	  ///
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern int testmode;

	  /// calling parameters
	  // extern Flags parameters;

	  // extern DLL_HEADER MeshingParameters mparam;

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern Array<int> tets_in_qualclass;

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern object tcl_todo_mutex;

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern volatile multithreadt multithread;

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern string ngdir;
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern DebugParameters debugparam;
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern bool verbose;

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern int h_argc;
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern char ** h_argv;


//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  DLL_HEADER extern weak_ptr<Mesh> global_mesh;
	  public static void SetGlobalMesh(Mesh m)
	  {
		PrintMessage(5, "set global mesh");
		global_mesh = m;
	  }

	  // stringstream emptystr;
	  // ostream * testout = &emptystr;
	  // testout -> clear(ios::failbit);

	  // ostream * testout = &cout;
	  public static ostream testout = new ostream(0);

	  // NetgenOutStream * testout = new NetgenOutStream;

	  public static ostream mycout = cout;
	  public static ostream myerr = cerr;

	  // some functions (visualization) still need a global mesh
	  // TraceGlobal glob1("global1");
	  public static DLL_HEADER Mesh * mesh = new DLL_HEADER();
	  public static DLL_HEADER NetgenGeometry * ng_geometry = new DLL_HEADER();
	  // TraceGlobal glob2("global2");

	  // global communicator for netgen
	  // DLL_HEADER NgMPI_Comm ng_comm;

	  public static weak_ptr<Mesh> global_mesh = new weak_ptr<Mesh>();

	  // true if netgen was started using the netgen executable
	  // false if netgen.gui was imported from python
	  public static DLL_HEADER bool netgen_executable_started = false;

	  //  Flags parameters;
	  public static int silentflag = 0;
	  public static int testmode = 0;

//C++ TO C# CONVERTER TODO TASK: 'volatile' has a different meaning in C#:
//ORIGINAL LINE: volatile multithreadt multithread;
	  public static multithreadt multithread = new multithreadt();

	  public static string ngdir = ".";

	  public static void Ng_PrintDest(string s)
	  {
		if (id == 0)
		{
		  mycout << s << flush;
		}
	  }

	  public static DLL_HEADER void MyError(string ch)
	  {
		Console.Write(ch);
		testout << "Error !!! " << ch << "\n" << flush;
	  }

	  internal static clock_t starttimea = new clock_t();



	  public static Array<int> tets_in_qualclass = new Array<int>();

	  public static object tcl_todo_mutex = new object();

	  public static int h_argc = 0;
	  public static string[] h_argv = null;

	  public static DebugParameters debugparam = new DebugParameters();
	  public static bool verbose = false;

	  public static uint timestamp = 0;
	public static string[] hexrules = {"rule \"Hexa left-right-top\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "flags t;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(1, 1, 0) { 1 };\n", \ "(0, 1, 0) { 1 };\n", \ "(0, 0, 1) { 1 };\n", \ "(1, 0, 1) { 1 };\n", \ "(1, 1, 1) { 1 };\n", \ "(0, 1, 1) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(4, 3, 2, 1) del;\n", \ "(3, 7, 6, 2) del;\n", \ "(7, 8, 5, 6) del;\n", \ "(8, 4, 1, 5) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 6, 2, 1);\n", \ "(7, 8, 4, 3);\n", \ "\n", \ "elements\n", \ "(4, 3, 2, 1, 8, 7, 6, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 0.3 P1, 0.3 P2, 0.3 P5, 0.3 P6, -0.05 P3, -0.05 P4, -0.05 P7, -0.05 P8 };\n", \ "{ 0.3 P3, 0.3 P4, 0.3 P7, 0.3 P8, -0.05 P1, -0.05 P2, -0.05 P5, -0.05 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 0.25 P1, 0.25 P2, 0.25 P5, 0.25 P6, -0.0 P3, -0.0 P4, -0.0 P7, -0.0 P8 };\n", \ "{ 0.25 P3, 0.25 P4, 0.25 P7, 0.25 P8, -0.0 P1, -0.0 P1, -0.0 P5, -0.0 P6 };\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Hexa left-right-top (10)\"\n", \ "\n", \ "quality 10\n", \ "\n", \ "flags t;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(1, 1, 0) { 1 };\n", \ "(0, 1, 0) { 1 };\n", \ "(0, 0, 1) { 1 };\n", \ "(1, 0, 1) { 1 };\n", \ "(1, 1, 1) { 1 };\n", \ "(0, 1, 1) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(4, 3, 2, 1) del;\n", \ "(3, 7, 6, 2) del;\n", \ "(7, 8, 5, 6) del;\n", \ "(8, 4, 1, 5) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 6, 2, 1);\n", \ "(7, 8, 4, 3);\n", \ "\n", \ "elements\n", \ "(4, 3, 2, 1, 8, 7, 6, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 0.251 P1, 0.251 P2, 0.251 P5, 0.251 P6, -0.05 P3, -0.001 P4, -0.001 P7, -0.001 P8 };\n", \ "{ 0.251 P3, 0.251 P4, 0.251 P7, 0.251 P8, -0.05 P1, -0.001 P2, -0.001 P5, -0.001 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 0.25 P1, 0.25 P2, 0.25 P5, 0.25 P6, -0.0 P3, -0.0 P4, -0.0 P7, -0.0 P8 };\n", \ "{ 0.25 P3, 0.25 P4, 0.25 P7, 0.25 P8, -0.0 P1, -0.0 P1, -0.0 P5, -0.0 P6 };\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Hexa left-right-top-front\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "flags t;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(1, 1, 0) { 1 };\n", \ "(0, 1, 0) { 1 };\n", \ "(0, 0, 1) { 1 };\n", \ "(1, 0, 1) { 1 };\n", \ "(1, 1, 1) { 1 };\n", \ "(0, 1, 1) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(4, 3, 2, 1) del;\n", \ "(3, 7, 6, 2) del;\n", \ "(7, 8, 5, 6) del;\n", \ "(8, 4, 1, 5) del;\n", \ "(1, 2, 6, 5) del;\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(7, 8, 4, 3);\n", \ "\n", \ "elements\n", \ "(4, 3, 2, 1, 8, 7, 6, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 0.3 P3, 0.3 P4, 0.3 P7, 0.3 P8, -0.05 P1, -0.05 P2, -0.05 P5, -0.05 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 0.25 P3, 0.25 P4, 0.25 P7, 0.25 P8, -0.0 P1, -0.0 P1, -0.0 P5, -0.0 P6 };\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Hexa fill\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "flags t;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(1, 1, 0) { 1 };\n", \ "(0, 1, 0) { 1 };\n", \ "(0, 0, 1) { 1 };\n", \ "(1, 0, 1) { 1 };\n", \ "(1, 1, 1) { 1 };\n", \ "(0, 1, 1) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(4, 3, 2, 1) del;\n", \ "(3, 7, 6, 2) del;\n", \ "(7, 8, 5, 6) del;\n", \ "(8, 4, 1, 5) del;\n", \ "(1, 2, 6, 5) del;\n", \ "(3, 4, 8, 7) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "\n", \ "elements\n", \ "(4, 3, 2, 1, 8, 7, 6, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P5 };\n", \ "{ 1 P3 };\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ 0};


	  public static HPRef_Struct Get_HPRef_Struct(HPREF_ELEMENT_TYPE type)
	  {
		HPRef_Struct hps = null;

		switch (type)
		{
		  case HPREF_ELEMENT_TYPE.HP_SEGM:
		hps = refsegm;
		break;
		  case HPREF_ELEMENT_TYPE.HP_SEGM_SINGCORNERL:
		hps = refsegm_scl;
		break;
		  case HPREF_ELEMENT_TYPE.HP_SEGM_SINGCORNERR:
		hps = refsegm_scr;
		break;
		  case HPREF_ELEMENT_TYPE.HP_SEGM_SINGCORNERS:
		hps = refsegm_sc2;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TRIG:
		hps = reftrig;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGCORNER:
		hps = reftrig_singcorner;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGCORNER12:
		hps = reftrig_singcorner12;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGCORNER123:
		hps = reftrig_singcorner123;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGCORNER123_2D:
		hps = reftrig_singcorner123_2D;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGE:
		hps = reftrig_singedge;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER1:
		hps = reftrig_singedgecorner1;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER2:
		hps = reftrig_singedgecorner2;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER12:
		hps = reftrig_singedgecorner12;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER3:
		hps = reftrig_singedgecorner3;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER13:
		hps = reftrig_singedgecorner13;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER23:
		hps = reftrig_singedgecorner23;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGECORNER123:
		hps = reftrig_singedgecorner123;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGES:
		hps = reftrig_singedges;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGES2:
		hps = reftrig_singedges2;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGES3:
		hps = reftrig_singedges3;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_SINGEDGES23:
		hps = reftrig_singedges23;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG_3SINGEDGES:
		hps = reftrig_3singedges;
		break;


		  case HPREF_ELEMENT_TYPE.HP_QUAD:
		hps = refquad;
		break;
		  case HPREF_ELEMENT_TYPE.HP_DUMMY_QUAD_SINGCORNER:
		hps = refdummyquad_singcorner;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_SINGCORNER:
		hps = refquad_singcorner;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_SINGEDGE:
		hps = refquad_singedge;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_0E_2VA:
		hps = refquad_0e_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_0E_2VB:
		hps = refquad_0e_2vb;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_0E_3V:
		hps = refquad_0e_3v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_0E_4V:
		hps = refquad_0e_4v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_1VA:
		hps = refquad_1e_1va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_1VB:
		hps = refquad_1e_1vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_1VC:
		hps = refquad_1e_1vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_1VD:
		hps = refquad_1e_1vd;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_2VA:
		hps = refquad_1e_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_2VB:
		hps = refquad_1e_2vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_2VC:
		hps = refquad_1e_2vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_2VD:
		hps = refquad_1e_2vd;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_2VE:
		hps = refquad_1e_2ve;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_2VF:
		hps = refquad_1e_2vf;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_3VA:
		hps = refquad_1e_3va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_3VB:
		hps = refquad_1e_3vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_3VC:
		hps = refquad_1e_3vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_3VD:
		hps = refquad_1e_3vd;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_1E_4V:
		hps = refquad_1e_4v;
		break;


		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E:
		hps = refquad_2e;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_1VA:
		hps = refquad_2e_1va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_1VB:
		hps = refquad_2e_1vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_1VC:
		hps = refquad_2e_1vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_2VA:
		hps = refquad_2e_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_2VB:
		hps = refquad_2e_2vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_2VC:
		hps = refquad_2e_2vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2E_3V:
		hps = refquad_2e_3v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_0V:
		hps = refquad_2eb_0v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_1VA:
		hps = refquad_2eb_1va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_1VB:
		hps = refquad_2eb_1vb;
		break;


		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_2VA:
		hps = refquad_2eb_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_2VB:
		hps = refquad_2eb_2vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_2VC:
		hps = refquad_2eb_2vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_2VD:
		hps = refquad_2eb_2vd;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_3VA:
		hps = refquad_2eb_3va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_3VB:
		hps = refquad_2eb_3vb;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_2EB_4V:
		hps = refquad_2eb_4v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_QUAD_3E:
		hps = refquad_3e;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_3E_3VA:
		hps = refquad_3e_3va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_3E_3VB:
		hps = refquad_3e_3vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD_3E_4V:
		hps = refquad_3e_4v;
		break;


		  case HPREF_ELEMENT_TYPE.HP_QUAD_4E:
		hps = refquad_4e;
		break;


		  case HPREF_ELEMENT_TYPE.HP_TET:
		hps = reftet;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_0E_1V:
		hps = reftet_0e_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_0E_2V:
		hps = reftet_0e_2v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_0E_3V:
		hps = reftet_0e_3v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_0E_4V:
		hps = reftet_0e_4v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TET_1E_0V:
		hps = reftet_1e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_1VA:
		hps = reftet_1e_1va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_1VB:
		hps = reftet_1e_1vb;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TET_1E_2VA:
		hps = reftet_1e_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_2VB:
		hps = reftet_1e_2vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_2VC:
		hps = reftet_1e_2vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_2VD:
		hps = reftet_1e_2vd;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TET_1E_3VA:
		hps = reftet_1e_3va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_3VB:
		hps = reftet_1e_3vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_4V:
		hps = reftet_1e_4v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_0V:
		hps = reftet_2ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_1VB:
		hps = reftet_2ea_1vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_1VC:
		hps = reftet_2ea_1vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_1VA:
		hps = reftet_2ea_1va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_2VA:
		hps = reftet_2ea_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_2VB:
		hps = reftet_2ea_2vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_2VC:
		hps = reftet_2ea_2vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EA_3V:
		hps = reftet_2ea_3v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_0V:
		hps = reftet_2eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_1V:
		hps = reftet_2eb_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_2VA:
		hps = reftet_2eb_2va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_2VB:
		hps = reftet_2eb_2vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_2VC:
		hps = reftet_2eb_2vc;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_3V:
		hps = reftet_2eb_3v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2EB_4V:
		hps = reftet_2eb_4v;
		break;


		  case HPREF_ELEMENT_TYPE.HP_TET_3EA_0V:
		hps = reftet_3ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EA_1V:
		hps = reftet_3ea_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EA_2V:
		hps = reftet_3ea_2v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EA_3V:
		hps = reftet_3ea_3v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_TET_3EB_0V:
		hps = reftet_3eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EB_1V:
		hps = reftet_3eb_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EB_2V:
		hps = reftet_3eb_2v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EC_0V:
		hps = reftet_3ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EC_1V:
		hps = reftet_3ec_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_3EC_2V:
		hps = reftet_3ec_2v;
		break;


		  case HPREF_ELEMENT_TYPE.HP_TET_1F_0E_0V:
		hps = reftet_1f_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1F_0E_1VA:
		hps = reftet_1f_0e_1va;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1F_0E_1VB:
		hps = reftet_1f_0e_1vb;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1F_1EA_0V:
		hps = reftet_1f_1ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_1F_1EB_0V:
		hps = reftet_1f_1eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_TET_2F_0E_0V:
		hps = reftet_2f_0e_0v;
		break;


		  case HPREF_ELEMENT_TYPE.HP_PRISM:
		hps = refprism;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_SINGEDGE:
		hps = refprism_singedge;
		break;
		//      case HP_PRISM_SINGEDGE_H1:
		//	hps = &refprism_singedge_h1; break;
		// case HP_PRISM_SINGEDGE_H12:
		//	hps = &refprism_singedge_h12; break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_SINGEDGE_V12:
		hps = refprism_singedge_v12;
		break;


		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_0E_0V:
		hps = refprism_1fa_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_0E_0V:
		hps = refprism_2fa_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FB_0E_0V:
		hps = refprism_1fb_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FB_1EA_0V:
		hps = refprism_1fb_1ea_0v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1E_0V:
		hps = refprism_1fa_1e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1E_0V:
		hps = refprism_2fa_1e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_0E_0V:
		hps = refprism_1fa_1fb_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1FB_0E_0V:
		hps = refprism_2fa_1fb_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_1EA_0V:
		hps = refprism_1fa_1fb_1ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_1EB_0V:
		hps = refprism_1fa_1fb_1eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1FB_1EA_0V:
		hps = refprism_2fa_1fb_1ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FB_1EC_0V:
		hps = refprism_1fb_1ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_1EC_0V:
		hps = refprism_1fa_1fb_1ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1FB_1EC_0V:
		hps = refprism_2fa_1fb_1ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FB_2EA_0V:
		hps = refprism_1fb_2ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_2EA_0V:
		hps = refprism_1fa_1fb_2ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1FB_2EA_0V:
		hps = refprism_2fa_1fb_2ea_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FB_2EB_0V:
		hps = refprism_1fb_2eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_2EB_0V:
		hps = refprism_1fa_1fb_2eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_2EC_0V:
		hps = refprism_1fa_1fb_2ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1FB_2EB_0V:
		hps = refprism_2fa_1fb_2eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FB_3E_0V:
		hps = refprism_1fb_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_1FB_3E_0V:
		hps = refprism_1fa_1fb_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_1FB_3E_0V:
			hps = refprism_2fa_1fb_3e_0v;
			break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FB_0E_0V:
		hps = refprism_2fb_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_2FB_0E_0V:
		hps = refprism_1fa_2fb_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_2FB_0E_0V:
			hps = refprism_2fa_2fb_0e_0v;
			break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FB_1EC_0V:
		hps = refprism_2fb_1ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_2FB_1EC_0V:
			hps = refprism_1fa_2fb_1ec_0v;
			break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_2FB_1EC_0V:
		hps = refprism_2fa_2fb_1ec_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_2FB_1EB_0V:
		hps = refprism_1fa_2fb_1eb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FB_3E_0V:
		hps = refprism_2fb_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_2FB_3E_0V:
		hps = refprism_1fa_2fb_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_2FB_3E_0V:
		hps = refprism_2fa_2fb_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_2E_0V:
		hps = refprism_1fa_2e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_2E_0V:
		hps = refprism_2fa_2e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_3E_0V:
		hps = refprism_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_3E_0V:
		hps = refprism_1fa_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_3E_0V:
		hps = refprism_2fa_3e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_3FB_0V:
		hps = refprism_3fb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_1FA_3FB_0V:
		hps = refprism_1fa_3fb_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM_2FA_3FB_0V:
		hps = refprism_2fa_3fb_0v;
		break;
		//  case HP_PRISM_3E_4EH:
		//  hps = &refprism_3e_4eh; break;   


		/*case HP_PRISM_1FB_1EB_0V:
		hps = &refprism_1fb_1eb_0v; break;
		  case HP_PRISM_2F_0E_0V:
		hps = &refprism_2f_0e_0v; break;
		*/


		  case HPREF_ELEMENT_TYPE.HP_PYRAMID:
		hps = refpyramid;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PYRAMID_0E_1V:
		hps = refpyramid_0e_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PYRAMID_EDGES:
		hps = refpyramid_edges;
		break;
		  case HPREF_ELEMENT_TYPE.HP_PYRAMID_1FB_0E_1VA:
		hps = refpyramid_1fb_0e_1va;
		break;


		  case HPREF_ELEMENT_TYPE.HP_HEX:
		hps = refhex;
		break;
		  case HPREF_ELEMENT_TYPE.HP_HEX_0E_1V:
		hps = refhex_0e_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_HEX_1E_1V:
		hps = refhex_1e_1v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_HEX_1E_0V:
		hps = refhex_1e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_HEX_3E_0V:
		hps = refhex_3e_0v;
		break;

		  case HPREF_ELEMENT_TYPE.HP_HEX_1F_0E_0V:
		hps = refhex_1f_0e_0v;
		break;
		  case HPREF_ELEMENT_TYPE.HP_HEX_1FA_1FB_0E_0V:
		hps = refhex_1fa_1fb_0e_0v;
		break;



		  default:
		  {
			  hps = null;
		  }
		break;
		}

		/*
		if (type != HP_TET_1E_4V && type != HP_TET_1E_2VD)
		  {
		if (hps->geom == HP_TET)
		  hps = &reftet;
		if (hps->geom == HP_TRIG)
		  hps = &reftrig;
		  }
		*/

		if (hps == null)
		{
		Console.Write("Attention hps : hp-refinement not implemented for case ");
		Console.Write(type);
		Console.Write("\n");
		PrintSysError("hp-refinement not implemented for case ", type);
		}

		return hps;
	  }

	  public static bool CheckSingularities(Mesh mesh, INDEX_2_HASHTABLE<int> edges, INDEX_2_HASHTABLE<int> edgepoint_dom, BitArray cornerpoint, BitArray edgepoint, INDEX_3_HASHTABLE<int> faces, INDEX_2_HASHTABLE<int> face_edges, INDEX_2_HASHTABLE<int> surf_edges, ref Array<int, PointIndex.BASE> facepoint, ref int levels, ref int act_ref)
	  {
		bool sing = false;
		if (mesh.GetDimension() == 3)
		{
		  /*
		  // check, if point has as least 3 different surfs:
	  
		  Array<INDEX_3, PointIndex::BASE> surfonpoint(mesh.GetNP());
			surfonpoint = INDEX_3(0,0,0);
	  
		  for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
		    {
		      const Element2d & el = mesh[sei];
		      int ind = el.GetIndex();
		      for (int j = 0; j < el.GetNP(); j++)
		        {
		  	INDEX_3 & i3 = surfonpoint[el[j]];
		  	if (ind != i3.I1() && ind != i3.I2() && ind != i3.I3())
		  	  {
		  	    i3.I1() = i3.I2();
		  	    i3.I2() = i3.I3();
		  	    i3.I3() = ind;
		  	  }
		        }
		    }
		  for (int i = 1; i <= mesh.GetNP(); i++)
		    if (surfonpoint.Get(i).I1())
		      cornerpoint.Set(i);
		  */
		  cornerpoint.Clear();

		  for (int i = 1; i <= mesh.GetNP(); i++)
		  {
			  if (new mesh.Point(i).Singularity() * levels >= act_ref != null)
			  {
			  cornerpoint.Set(i);
			  sing = true;
			  }
		  }
		  Console.Write("\n");

		  for (int i = 1; i <= mesh.GetNSeg(); i++)
		  {
			if (mesh.LineSegment(i).singedge_left * levels >= act_ref)
			{
				INDEX_2 i2 = new INDEX_2(mesh.LineSegment i[0], mesh.LineSegment i[1]);

				/*
			  // before
				edges.Set (i2, 1);
				i2.Sort();
				INDEX_2 i2s(i2.I2(), i2.I1());
				edges.Set (i2s, 1);
				*/

				edges.Set(i2, 1);
				INDEX_2 i2s = new INDEX_2(i2.I2(), i2.I1());
				edges.Set(i2s, 1);


				edgepoint.Set(i2.I1());
				edgepoint.Set(i2.I2());
				sing = true;
			}
		  }

		  // if 2 adjacent edges of an element are singular, the
		  // commen point must be a singular point
		  for (int i = 1; i <= mesh.GetNE(); i++)
		  {
			  Element el = mesh.VolumeElement(i);
			  ELEMENT_EDGE[] eledges = MeshTopology.GetEdges1(el.GetType());
			  int nedges = MeshTopology.GetNEdges(el.GetType());
			  for (int j = 0; j < nedges; j++)
			  {
				for (int k = 0; k < nedges; k++)
				{
			  if (j != k)
			  {
				  INDEX_2 ej = new INDEX_2(el.PNum(eledges[j][0]), el.PNum(eledges[j][1]));
				  ej.Sort();
				  INDEX_2 ek = new INDEX_2(el.PNum(eledges[k][0]), el.PNum(eledges[k][1]));
				  ek.Sort();
				  if (edges.Used(ej) && edges.Used(ek))
				  {
				  if (ej.I1() == ek.I1())
				  {
					  cornerpoint.Set(ek.I1());
				  }
				  if (ej.I1() == ek.I2())
				  {
					  cornerpoint.Set(ek.I2());
				  }
				  if (ej.I2() == ek.I1())
				  {
					  cornerpoint.Set(ek.I1());
				  }
				  if (ej.I2() == ek.I2())
				  {
					  cornerpoint.Set(ek.I2());
				  }
				  }
			  }
				}
			  }
		  }

		  edgepoint.Or(cornerpoint);
		  (*testout) << "cornerpoint = " << "\n" << cornerpoint << "\n";
		  (*testout) << "edgepoint = " << "\n" << edgepoint << "\n";

		  facepoint = 0;
		  for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
		  {
			  Element2d el = mesh[sei];
			  FaceDescriptor fd = mesh.GetFaceDescriptor(el.GetIndex());

			  int domnr = 0;
			  if (fd.DomainInSingular() * levels < act_ref && fd.DomainOutSingular() * levels < act_ref)
			  {
					domnr = 0;
					continue;
			  }

			  if (fd.DomainInSingular() * levels >= act_ref)
			  {
			  domnr = fd.DomainIn();
			  sing = true;
			  }
			  if (fd.DomainOutSingular() * levels >= act_ref)
			  {
			  domnr = fd.DomainOut();
			  sing = true;
			  }
			  if (fd.DomainInSingular() * levels >= act_ref && fd.DomainOutSingular() * levels >= act_ref)
			  {
			  domnr = -1;
			  sing = true;
			  }

			  INDEX_3 i3 = new INDEX_3();
			  if (el.GetNP() == 3)
			  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: i3 = INDEX_3::Sort(el[0], el[1], el[2]);
				i3.CopyFrom(INDEX_3.Sort(new netgen.Element2d(el[0]), new netgen.Element2d(el[1]), new netgen.Element2d(el[2])));
			  }
			  else
			  {
			  INDEX_4 i4 = new INDEX_4(el[0], el[1], el[2], el[3]);
			  i4.Sort();
			  i3 = new INDEX_3(i4.I1(), i4.I2(), i4.I3());
			  }
			  faces.Set(i3, domnr);

			  for (int j = 0; j < el.GetNP(); j++)
			  {
			  face_edges.Set(INDEX_2.Sort(new netgen.Element2d(el[j]), new netgen.Element2d(el[(j + 1) % el.GetNP()])), domnr);

			  surf_edges.Set(INDEX_2.Sort(new netgen.Element2d(el[j]), new netgen.Element2d(el[(j + 1) % el.GetNP()])), fd.SurfNr() + 1);

			  facepoint[el[j]] = domnr;
			  }

		  }
		  (*testout) << "singular faces = " << faces << "\n";
		  (*testout) << "singular faces_edges = " << face_edges << "\n";
		}
		  else
		  {
		  // 2D case

		  // check, if point has as least 3 different surfs:
		  Array<INDEX_3, PointIndex.BASE> surfonpoint = new Array<INDEX_3, PointIndex.BASE>(mesh.GetNP());

		  for (int i = 1; i <= mesh.GetNP(); i++)
		  {
			surfonpoint.Elem(i) = new INDEX_3(0, 0, 0);
		  }

		  for (int i = 1; i <= mesh.GetNSeg(); i++)
		  {
			  Segment seg = mesh.LineSegment(i);
			  int ind = seg.edgenr;

			  if (seg.singedge_left * levels >= act_ref)
			  {
			  INDEX_2 i2 = new INDEX_2(mesh.LineSegment i[0], mesh.LineSegment i[1]);
			  edges.Set(i2, 1);
			  edgepoint.Set(i2.I1());
			  edgepoint.Set(i2.I2());
			  *testout << " singleft " << "\n";
			  *testout << " mesh.LineSegment(i).domout " << mesh.LineSegment(i).domout << "\n";
			  *testout << " mesh.LineSegment(i).domin " << mesh.LineSegment(i).domin << "\n";
			  edgepoint_dom.Set(new INDEX_2(mesh.LineSegment(i).domin, i2.I1()), 1);
			  edgepoint_dom.Set(new INDEX_2(mesh.LineSegment(i).domin, i2.I2()), 1);
			  sing = true;

			  }

			  if (seg.singedge_right * levels >= act_ref)
			  {
			  INDEX_2 i2 = new INDEX_2(mesh.LineSegment i[1], mesh.LineSegment i[0]);
			  edges.Set(i2, 1);
			  edgepoint.Set(i2.I1());
			  edgepoint.Set(i2.I2());

			  *testout << " singright " << "\n";
			  *testout << " mesh.LineSegment(i).domout " << mesh.LineSegment(i).domout << "\n";
			  *testout << " mesh.LineSegment(i).domin " << mesh.LineSegment(i).domin << "\n";

			  edgepoint_dom.Set(new INDEX_2(mesh.LineSegment(i).domout, i2.I1()), 1);
			  edgepoint_dom.Set(new INDEX_2(mesh.LineSegment(i).domout, i2.I2()), 1);
			  sing = true;
			  }

			  // (*testout) << "seg = " << ind << ", " << seg[0] << "-" << seg[1] << endl;


			  if (seg.singedge_left * levels >= act_ref || seg.singedge_right * levels >= act_ref)
			  {
			  for (int j = 0; j < 2; j++)
			  {
				  int pi = (j == 0) ? seg[0] : seg[1];
				  INDEX_3 i3 = surfonpoint.Elem(pi);
				  if (ind != i3.I1() && ind != i3.I2())
				  {
				  i3.I1() = i3.I2();
				  i3.I2() = ind;
				  }
			  }
			  }
		  }


		  for (int i = 1; i <= mesh.GetNP(); i++)
		  {
			  // mark points for refinement that are in corners between two anisotropic edges
			  if (surfonpoint.Get(i).I1())
			  {
			  // cornerpoint.Set(i);    // disabled by JS, Aug 2009
			  edgepoint.Set(i);
			  }

			  // mark points for refinement that are explicitly specified in input file
			  if (new mesh.Point(i).Singularity() * levels >= act_ref != null)
			  {
			  cornerpoint.Set(i);
			  edgepoint.Set(i);
			  sing = true;
			  }
		  }

		  edgepoint.Or(cornerpoint);

		  (*testout) << "2d sing edges: " << "\n" << edges << "\n";
		  (*testout) << "2d cornerpoints: " << "\n" << cornerpoint << "\n" << "2d edgepoints: " << "\n" << edgepoint << "\n";

		  facepoint = 0;
		  }

		  if (!sing)
		  {
			Console.Write("PrepareElements no more to do for actual refinement ");
			Console.Write(act_ref);
			Console.Write("\n");
		  }

		  return (sing);
	  }

	  public static bool ClassifyHPElements(Mesh mesh, Array<HPRefElement> elements, ref int act_ref, ref int levels)
	  {
		INDEX_2_HASHTABLE<int> edges = new INDEX_2_HASHTABLE<int>(mesh.GetNSeg() + 1);
		BitArray edgepoint = new BitArray(mesh.GetNP());
		INDEX_2_HASHTABLE<int> edgepoint_dom = new INDEX_2_HASHTABLE<int>(mesh.GetNSeg() + 1);

		edgepoint.Clear();
		BitArray cornerpoint = new BitArray(mesh.GetNP());
		cornerpoint.Clear();

		// value = nr > 0 ... refine elements in domain nr
		// value = -1   ..... refine elements in any domain
		INDEX_3_HASHTABLE<int> faces = new INDEX_3_HASHTABLE<int>(mesh.GetNSE() + 1);
		INDEX_2_HASHTABLE<int> face_edges = new INDEX_2_HASHTABLE<int>(mesh.GetNSE() + 1);
		INDEX_2_HASHTABLE<int> surf_edges = new INDEX_2_HASHTABLE<int>(mesh.GetNSE() + 1);
		Array<int, PointIndex.BASE> facepoint = new Array<int, PointIndex.BASE>(mesh.GetNP());

		bool sing = CheckSingularities(mesh, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, ref facepoint, ref levels, ref act_ref);

		if (sing == false)
		{
			return (sing);
		}

		int cnt_undef = 0;
		int cnt_nonimplement = 0;
		Array<int> misses = new Array<int>(10000);
		misses = 0;

		(*testout) << "edgepoint_dom = " << "\n" << edgepoint_dom << "\n";

		for (int i = 0; i < elements.Size(); i++)
		{
		// *testout << "classify element " << i << endl;

		HPRefElement hpel = elements[i];
		HPRef_Struct hprs = Get_HPRef_Struct(hpel.type);
		HPRefElement old_el = elements[i];
		int dd = 3;


		if (act_ref != 1 && (hpel.type == HPREF_ELEMENT_TYPE.HP_HEX || hpel.type == HPREF_ELEMENT_TYPE.HP_PRISM || hpel.type == HPREF_ELEMENT_TYPE.HP_TET || hpel.type == HPREF_ELEMENT_TYPE.HP_PYRAMID || hpel.type == HPREF_ELEMENT_TYPE.HP_QUAD || hpel.type == HPREF_ELEMENT_TYPE.HP_TRIG || hpel.type == HPREF_ELEMENT_TYPE.HP_SEGM))
		{
		  continue;
		}

		sing = true;
		switch (hprs.geom)
		{
		  case HPREF_ELEMENT_TYPE.HP_TET:
		  {
			  hpel.type = ClassifyTet(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint);
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_PRISM:
		  {
			  hpel.type = ClassifyPrism(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint);


			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_HEX:
		  {
			  hpel.type = ClassifyHex(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint);
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_TRIG:
		  {
			  int dim = mesh.GetDimension();
			  FaceDescriptor fd = mesh.GetFaceDescriptor(hpel.GetIndex());

			  hpel.type = ClassifyTrig(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint, dim, fd);

			  dd = 2;
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_QUAD:
		  {
			  int dim = mesh.GetDimension();
			  FaceDescriptor fd = mesh.GetFaceDescriptor(hpel.GetIndex());
			  hpel.type = ClassifyQuad(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint, dim, fd);

			  dd = 2;
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_SEGM:
		  {
			  hpel.type = ClassifySegm(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint);
			  dd = 1;
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_PYRAMID:
		  {
			  hpel.type = ClassifyPyramid(hpel, edges, edgepoint_dom, cornerpoint, edgepoint, faces, face_edges, surf_edges, facepoint);

			  Console.Write(" ** Pyramid classified  ");
			  Console.Write(hpel.type);
			  Console.Write("\n");
			  break;
		  }
		  default:
		  {
			  Console.Write("illegal element type for hp-prepare elements ");
			  Console.Write(hpel.type);
			  Console.Write("\n");
			  throw new Exception("hprefinement.cpp: don't know how to set parameters");
		  }
		}

		if (hpel.type == HPREF_ELEMENT_TYPE.HP_NONE)
		{
		  cnt_undef++;
		}

		//else
		//cout << "elem " << i << " classified type " << hpel.type << endl;



		if (!Get_HPRef_Struct(hpel.type))
		{
			(*testout) << "hp-element-type " << hpel.type << " not implemented   " << "\n";
			(*testout) << " elType " << hprs.geom << "\n";
	 (cout) << " elType " << hprs.geom << "\n";
			cnt_nonimplement++;
			misses[(int)hpel.type]++;
		}


		for (int j = 0; j < hpel.np; j++)
		{
			for (int k = 0; k < hpel.np; k++)
			{
			  if (hpel[j] == old_el.pnums[k])
			  {
			  for (int l = 0;l < dd;l++)
			  {
				hpel.param[j][l] = old_el.param[k][l];
			  }
			  break;
			  }
			}
		}

		}


		Console.Write("undefined elements update classification: ");
		Console.Write(cnt_undef);
		Console.Write("\n");
		Console.Write("non-implemented in update classification: ");
		Console.Write(cnt_nonimplement);
		Console.Write("\n");

		for (int i = 0; i < misses.Size(); i++)
		{
		  if (misses[i])
		  {
		Console.Write(" in update classification missing case ");
		Console.Write(i);
		Console.Write(" occurred ");
		Console.Write(misses[i]);
		Console.Write(" times");
		Console.Write("\n");
		  }
		}

		return (sing);
	  }


	  public static void InitHPElements(Mesh mesh, Array<HPRefElement> elements)
	  {
		for (ElementIndex i = 0; i < mesh.GetNE(); i++)
		{
		HPRefElement hpel = new HPRefElement(mesh[i]);
		hpel.coarse_elnr = i;

		switch (mesh[i].GetType())
		{
		  case ELEMENT_TYPE.PRISM:
			  hpel.type = HPREF_ELEMENT_TYPE.HP_PRISM;
			  break;
		  case ELEMENT_TYPE.HEX:
			  hpel.type = HPREF_ELEMENT_TYPE.HP_HEX;
			  break;
		  case ELEMENT_TYPE.TET:
			  hpel.type = HPREF_ELEMENT_TYPE.HP_TET;
			  break;
		  case ELEMENT_TYPE.PYRAMID:
			  hpel.type = HPREF_ELEMENT_TYPE.HP_PYRAMID;
			  break;

			  default:
				cerr << "HPRefElement: illegal elementtype (1) " << mesh[i].GetType() << "\n";
				throw new Exception("HPRefElement: illegal elementtype (1)");
		}
		elements.Append(hpel);
		}

		for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
		{
		HPRefElement hpel = new HPRefElement(mesh[i]);
		hpel.coarse_elnr = i;
		switch (mesh[i].GetType())
		{
		  case ELEMENT_TYPE.TRIG:
			  hpel.type = HPREF_ELEMENT_TYPE.HP_TRIG;
			  break;
		  case ELEMENT_TYPE.QUAD:
			  hpel.type = HPREF_ELEMENT_TYPE.HP_QUAD;
			  break;

			  default:
				cerr << "HPRefElement: illegal elementtype (1b) " << mesh[i].GetType() << "\n";
				throw new Exception("HPRefElement: illegal elementtype (1b)");
		}
		elements.Append(hpel);
		}

		for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
		{
		Segment seg = mesh[i];
		HPRefElement hpel = new HPRefElement(mesh[i]);
		hpel.coarse_elnr = i;
		hpel.type = HPREF_ELEMENT_TYPE.HP_SEGM;
		hpel.index = seg.edgenr + 10000 * seg.si;
		if (seg.edgenr >= 10000)
		{
			throw new Exception("assumption that seg.edgenr < 10000 is wrong");
		}
		elements.Append(hpel);
		}
	  }



	  /* *******************************  DoRefinement *************************************** */
	  public static void DoRefinement(Mesh mesh, Array<HPRefElement> elements, Refinement @ref, double fac1)
	  {
		elements.SetAllocSize(5 * elements.Size());
		INDEX_2_HASHTABLE<int> newpts = new INDEX_2_HASHTABLE<int>(elements.Size() + 1);
		INDEX_3_HASHTABLE<int> newfacepts = new INDEX_3_HASHTABLE<int>(elements.Size() + 1);

		// prepare new points  

		fac1 = Math.Max(0.001,Math.Min(0.33,fac1));
		Console.Write(" in HP-REFINEMENT with fac1 ");
		Console.Write(fac1);
		Console.Write("\n");
		*testout << " in HP-REFINEMENT with fac1 " << fac1 << "\n";


		int oldelsize = elements.Size();

		for (int i = 0; i < oldelsize; i++)
		{
		HPRefElement el = elements[i];
		HPRef_Struct hprs = Get_HPRef_Struct(el.type);

		if (hprs == null)
		{
			Console.Write("Refinementstruct not defined for element ");
			Console.Write(el.type);
			Console.Write("\n");
			continue;
		}

		int j = 0;
		while (hprs.splitedges[j][0] != 0)
		{
			INDEX_2 i2 = new INDEX_2(el.pnums[hprs.splitedges[j][0] - 1], el.pnums[hprs.splitedges[j][1] - 1]);
			if (!newpts.Used(i2))
			{
			Point < 3> np;
			for (int l = 0;l < 3;l++)
			{
			  np(l) = (1 - fac1) * new mesh.Point(i2.I1())(l) + fac1 * new mesh.Point(i2.I2())(l);
			}

			int npi = mesh.AddPoint(np);
			newpts.Set(i2, npi);
			}
			j++;
		}

		j = 0;
		if (hprs.splitfaces)
		{
		  while (hprs.splitfaces[j][0] != 0)
		  {
			  INDEX_3 i3 = new INDEX_3(el.pnums[hprs.splitfaces[j][0] - 1], el.pnums[hprs.splitfaces[j][1] - 1], el.pnums[hprs.splitfaces[j][2] - 1]);

			  if (i3.I2() > i3.I3())
			  {
				  Swap(ref i3.I2(), ref i3.I3());
			  }

			  if (!newfacepts.Used(i3))
			  {
			  Point < 3> np;
				  for (int l = 0;l < 3;l++)
				  {
				  np(l) = (1 - 2 * fac1) * new mesh.Point(i3.I1())(l) + fac1 * new mesh.Point(i3.I2())(l) + fac1 * new mesh.Point(i3.I3())(l);
				  }
			  int npi = mesh.AddPoint(np);
			  newfacepts.Set(i3, npi);
			  }
			  j++;
		  }
		}
		}

		for (int i = 0; i < oldelsize; i++)
		{
		HPRefElement el = elements[i];
		HPRef_Struct hprs = Get_HPRef_Struct(el.type);
		int newlevel = el.levelx;
		int oldnp = 0;
		switch (hprs.geom)
		{
		  case HPREF_ELEMENT_TYPE.HP_SEGM:
			  oldnp = 2;
			  break;
		  case HPREF_ELEMENT_TYPE.HP_TRIG:
			  oldnp = 3;
			  break;
		  case HPREF_ELEMENT_TYPE.HP_QUAD:
			  oldnp = 4;
			  break;
		  case HPREF_ELEMENT_TYPE.HP_TET:
			  oldnp = 4;
			  break;
		  case HPREF_ELEMENT_TYPE.HP_PYRAMID:
			  oldnp = 5;
			  break;
		  case HPREF_ELEMENT_TYPE.HP_PRISM:
			  oldnp = 6;
			  break;
		  case HPREF_ELEMENT_TYPE.HP_HEX:
			  oldnp = 8;
			  break;

			  default:
				cerr << "HPRefElement: illegal type (3) " << hprs.geom << "\n";
				throw new Exception("HPRefElement::SetType: illegal type (3)");
		}


		if (el.type == HPREF_ELEMENT_TYPE.HP_SEGM || el.type == HPREF_ELEMENT_TYPE.HP_TRIG || el.type == HPREF_ELEMENT_TYPE.HP_QUAD || el.type == HPREF_ELEMENT_TYPE.HP_TET || el.type == HPREF_ELEMENT_TYPE.HP_PRISM || el.type == HPREF_ELEMENT_TYPE.HP_HEX || el.type == HPREF_ELEMENT_TYPE.HP_PYRAMID)
		{
		  newlevel = el.levelx;
		}

		if (hprs == null)
		{
			continue;
		}

		int[] newpnums = new int[64];
		double[][] newparam = RectangularArrays.RectangularDoubleArray(64, 3);

		int j;
		for (j = 0; j < oldnp; j++)
		{
			newpnums[j] = el.pnums[j];
			for (int l = 0; l < 3; l++)
			{
			  newparam[j][l] = el.param[j][l];
			}
		}

		// split edges, incl. transferring curvature
		j = 0;
		while (hprs.splitedges[j][0] != 0)
		{
			INDEX_2 i2 = new INDEX_2(el.pnums[hprs.splitedges[j][0] - 1], el.pnums[hprs.splitedges[j][1] - 1]);

			int npi = newpts.Get(i2);
			newpnums[hprs.splitedges[j][2] - 1] = npi;

			for (int l = 0; l < 3; l++)
			{
			  newparam[hprs.splitedges[j][2] - 1][l] = (1 - fac1) * el.param[hprs.splitedges[j][0] - 1][l] + fac1 * el.param[hprs.splitedges[j][1] - 1][l];
			}

			j++;
		}

		// split faces
		j = 0;
		if (hprs.splitfaces)
		{
		  while (hprs.splitfaces[j][0] != 0)
		  {
			  INDEX_3 i3 = new INDEX_3(el.pnums[hprs.splitfaces[j][0] - 1], el.pnums[hprs.splitfaces[j][1] - 1], el.pnums[hprs.splitfaces[j][2] - 1]);
			  if (i3.I2() > i3.I3())
			  {
			Swap(ref i3.I2(), ref i3.I3());
			  }
			  int npi = newfacepts.Get(i3);
			  newpnums[hprs.splitfaces[j][3] - 1] = npi;


			  for (int l = 0; l < 3; l++)
			  {
			newparam[hprs.splitfaces[j][3] - 1][l] = (1 - 2 * fac1) * el.param[hprs.splitfaces[j][0] - 1][l] + fac1 * el.param[hprs.splitfaces[j][1] - 1][l] + fac1 * el.param[hprs.splitfaces[j][2] - 1][l];
			  }
			  j++;
		  }
		}
		// split elements
		j = 0;
		if (hprs.splitelements)
		{
		  while (hprs.splitelements[j][0] != 0)
		  {
			  //int pi1 = el.pnums[hprs->splitelements[j][0]-1];
			  Point < 3> np;
				  for (int l = 0;l < 3;l++)
				  {
			  np(l) = (1 - 3 * fac1) * new mesh.Point(el.pnums[hprs.splitelements[j][0] - 1])(l) + fac1 * new mesh.Point(el.pnums[hprs.splitelements[j][1] - 1])(l) + fac1 * new mesh.Point(el.pnums[hprs.splitelements[j][2] - 1])(l) + fac1 * new mesh.Point(el.pnums[hprs.splitelements[j][3] - 1])(l);
				  }

			  int npi = mesh.AddPoint(np);

			  newpnums[hprs.splitelements[j][4] - 1] = npi;


			  for (int l = 0; l < 3; l++)
			  {
			newparam[hprs.splitelements[j][4] - 1][l] = (1 - 3 * fac1) * el.param[hprs.splitelements[j][0] - 1][l] + fac1 * el.param[hprs.splitelements[j][1] - 1][l] + fac1 * el.param[hprs.splitelements[j][2] - 1][l] + fac1 * el.param[hprs.splitelements[j][3] - 1][l];
			  }

			  j++;
		  }
		}

		j = 0;

		/*
		*testout << " newpnums = ";
		for (int hi = 0; hi < 64; hi++)
		  *testout << newpnums[hi] << " ";
		*testout << endl;
		*/

		while ((int)hprs.neweltypes[j] != 0)
		{
			HPRef_Struct hprsnew = Get_HPRef_Struct(hprs.neweltypes[j]);
			HPRefElement newel = new HPRefElement(el);

			newel.type = hprs.neweltypes[j];
			// newel.index = elements[i].index;
			// newel.coarse_elnr = elements[i].coarse_elnr;
				if (newel.type == HPREF_ELEMENT_TYPE.HP_SEGM || newel.type == HPREF_ELEMENT_TYPE.HP_TRIG || newel.type == HPREF_ELEMENT_TYPE.HP_QUAD || newel.type == HPREF_ELEMENT_TYPE.HP_TET || newel.type == HPREF_ELEMENT_TYPE.HP_PRISM || newel.type == HPREF_ELEMENT_TYPE.HP_HEX || newel.type == HPREF_ELEMENT_TYPE.HP_PYRAMID)
				{
				  newel.levelx = newel.levely = newel.levelz = newlevel;
				}
				else
				{
				  newel.levelx = newel.levely = newel.levelz = newlevel + 1;
				}

				switch (hprsnew.geom)
				{
			  case HPREF_ELEMENT_TYPE.HP_SEGM:
				  newel.np = 2;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_QUAD:
				  newel.np = 4;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_TRIG:
				  newel.np = 3;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_HEX:
				  newel.np = 8;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_PRISM:
				  newel.np = 6;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_TET:
				  newel.np = 4;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_PYRAMID:
				  newel.np = 5;
				  break;
				  default:
					throw new Exception("hprefinement.cpp: illegal type");
				}

			for (int k = 0; k < newel.np; k++)
			{
			  newel.pnums[k] = newpnums[hprs.newels[j][k] - 1];
			}

			/*
			*testout  << " newel pnums " ; 
			for (int k = 0; k < newel.np; k++)  
			  *testout  << newel.pnums[k] << "\t"; 
			*testout << endl; 
			*/

			for (int k = 0; k < newel.np; k++)
			{
			for (int l = 0; l < 3; l++)
			{
				newel.param[k][l] = newparam[hprs.newels[j][k] - 1][l];
				//    *testout << newel.param[k][l] << " \t ";
			}
			// *testout << endl; 
			}

			if (j == 0)
			{
			  elements[i] = newel; // overwrite old element
			}
			else
			{
			  elements.Append(newel);
			}
			j++;
		}
		}
	  }






	  /* ************************** DoRefineDummies ******************************** */

	  public static void DoRefineDummies(Mesh mesh, Array<HPRefElement> elements, Refinement @ref)
	  {
		int oldelsize = elements.Size();

		for (int i = 0; i < oldelsize; i++)
		{
		HPRefElement el = elements[i];

		HPRef_Struct hprs = Get_HPRef_Struct(el.type);
		if (hprs == null)
		{
			continue;
		}

		if (el.type != HPREF_ELEMENT_TYPE.HP_DUMMY_QUAD_SINGCORNER && el.type != HPREF_ELEMENT_TYPE.HP_PYRAMID_EDGES && el.type != HPREF_ELEMENT_TYPE.HP_PYRAMID_0E_1V && el.type != HPREF_ELEMENT_TYPE.HP_HEX_0E_1V && el.type != HPREF_ELEMENT_TYPE.HP_HEX_1E_1V && el.type != HPREF_ELEMENT_TYPE.HP_HEX_1E_0V && el.type != HPREF_ELEMENT_TYPE.HP_HEX_3E_0V)
		{
			continue;
		}

		int newlevel = el.levelx;

		int[] newpnums = new int[8];
		int j;
		for (j = 0; j < 8; j++)
		{
		  newpnums[j] = el.pnums[j];
		}

		double[][] newparam = RectangularArrays.RectangularDoubleArray(8, 3);
		for (j = 0; j < 8; j++)
		{
		  for (int k = 0; k < 3; k++)
		  {
			newparam[j][k] = el.param[j][k];
		  }
		}

		j = 0;
		while ((int)hprs.neweltypes[j] != 0)
		{
			HPRef_Struct hprsnew = Get_HPRef_Struct(hprs.neweltypes[j]);
			HPRefElement newel = new HPRefElement(el);
			switch (hprsnew.geom)
			{
			  case HPREF_ELEMENT_TYPE.HP_SEGM:
				  newel.np = 2;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_QUAD:
				  newel.np = 4;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_TRIG:
				  newel.np = 3;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_HEX:
				  newel.np = 8;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_PRISM:
				  newel.np = 6;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_TET:
				  newel.np = 4;
				  break;
			  case HPREF_ELEMENT_TYPE.HP_PYRAMID:
				  newel.np = 5;
				  break;

				  default:
					cerr << "HPRefElement: illegal type (4) " << hprsnew.geom << "\n";
					throw new Exception("HPRefElement: illegal type (4)");

			}
			newel.type = hprs.neweltypes[j];
			for (int k = 0; k < 8; k++)
			{
			  newel.pnums[k] = newpnums[hprs.newels[j][k] - 1];
			}
			newel.index = el.index;
			newel.coarse_elnr = el.coarse_elnr;
			newel.levelx = newel.levely = newel.levelz = newlevel;

			for (int k = 0; k < 8; k++)
			{
			  for (int l = 0; l < 3; l++)
			  {
			newel.param[k][l] = newparam[hprs.newels[j][k] - 1][l];
			  }
			}

			if (j == 0)
			{
			  elements[i] = newel;
			}
			else
			{
			  elements.Append(newel);
			}
			j++;
		}
		}
	  }







	  public static void SubdivideDegeneratedHexes(Mesh mesh, Array<HPRefElement> elements, double fac1)
	  {
		int oldne = elements.Size();
		for (int i = 0; i < oldne; i++)
		{
		  if (Get_HPRef_Struct(elements[i].type).geom == HPREF_ELEMENT_TYPE.HP_HEX)
		  {
		  bool common = false;
		  for (int j = 0; j < 8; j++)
		  {
			for (int k = 0; k < j; k++)
			{
			  if (elements[i].pnums[j] == elements[i].pnums[k])
			  {
			common = true;
			  }
			}
		  }
		  if (common)
		  {


			  Console.Write(" Degenerate Hex found ");
			  Console.Write("\n");
				  *testout << " Degenerate Hex found " << "\n";
			  HPRefElement el = elements[i];
			  HPRefElement newel = new HPRefElement(el);

			  Point < 3> center(0,0,0);
			  double[] newparam = {0, 0, 0};

			  for (int j = 0; j < 8; j++)
			  {


			  center += 0.125 * Vec < 3>(mesh[el.pnums[j]]);
			  // 0.125 originates form 8 points not from fac1;

			  for (int l = 0; l < 3; l++)
			  {
				newparam[l] += 0.125 * el.param[j][l];
			  }

			  }

			  int npi = mesh.AddPoint(center);

			  ELEMENT_FACE[] faces = MeshTopology.GetFaces1(ELEMENT_TYPE.HEX);

			  for (int j = 0; j < 6; j++)
			  {
			  Array<int> pts = new Array<int>();
			  for (int k = 0; k < 4; k++)
			  {
				  bool same = false;
				  for (int l = 0; l < pts.Size(); l++)
				  {
				if (el.pnums[pts[l]] == el.pnums[faces[j][k] - 1])
				{
				  same = true;
				}
				  }
				  if (!same)
				  {
				pts.Append(faces[j][k] - 1);
				  }

			  }


			  if (pts.Size() == 3) // TrigFace -> TET
			  {

				  for (int k = 0; k < 3; k++)
				  {
				  newel.pnums[k] = el.pnums[pts[2 - k]];
				  for (int l = 0; l < 3; l++)
				  {
					newel.param[k][l] = el.param[pts[2 - k]][l];
				  }
				  }
				  newel.pnums[3] = npi;
				  for (int l = 0; l < 3; l++)
				  {
				newel.param[3][l] = newparam[l];
				  }

				  newel.type = HPREF_ELEMENT_TYPE.HP_TET;
				  newel.np = 4;
			  }
			  else
			  {
				  for (int k = 0; k < 4; k++)
				  {
				  newel.pnums[k] = el.pnums[pts[3 - k]];
				  for (int l = 0; l < 3; l++)
				  {
					newel.param[k][l] = el.param[pts[3 - k]][l];
				  }
				  }

				  newel.pnums[4] = npi;
				  for (int l = 0; l < 3; l++)
				  {
				newel.param[4][l] = newparam[l];
				  }

				  newel.type = HPREF_ELEMENT_TYPE.HP_PYRAMID;
				  newel.np = 5;
			  }

			  if (j == 0)
			  {
				elements[i] = newel;
			  }
			  else
			  {
				elements.Append(newel);
			  }


			  }

			  /*     const ELEMENT_EDGE * edges = MeshTopology::GetEdges (HEX);
			   
			for(int k=0;k<12;k++) 
			  { 
				int e[2];  
				for(int l=0;l<2;l++) e[l] = edges[k][l]-1; 
				if(el.PNum(e[0]+1)!=el.PNum(e[1]+1)) 
				  { 
				newel.SetType(HP_SEGM);
				for(int l=0;l<2;l++) 
				  { 
					newel.pnums[0] = el.PNum(e[l]+1); 
					newel.pnums[1] = npi; 
					for(int j=0;j<3;j++) 
					  {
					//	newel.param[0][j] = el.param[e[l]][j]; 
					//	newel.param[1][j] = newparam[j]; 
					  } 
	
					elements.Append(newel);
				  }
				newel.SetType(HP_TRIG);
				newel.pnums[0] = el.PNum(e[0]+1); 			
				newel.pnums[1] = el.PNum(e[1]+1); 			
				newel.pnums[2] = npi; 
	
				*testout << "DEGHEX TRIG :: newpnums " << newel.pnums[0] << "\t"  << newel.pnums[1] << "\t"  << newel.pnums[2] << endl;  
		cout << "DEGHEX TRIG :: newpnums " << newel.pnums[0] << "\t"  << newel.pnums[1] << "\t"  << newel.pnums[2] << endl;  
				for(int j=0;j<3;j++) 
				  {
					// newel.param[0][j] = el.param[e[0]][j]; 
					//   newel.param[1][j] = el.param[e[1]][j]; 
					//   newel.param[2][j] = newparam[j]; 
				  } 
	
				elements.Append(newel);
				  }
	
				  }*/
		  }
		  }
		}
	  }


	  public static void CalcStatistics(Array<HPRefElement> elements)
	  {
		return;
	#if ABC
		int i;
		int p;
		int nsegm = 0;
		int ntrig = 0;
		int nquad = 0;
		int nhex = 0;
		int nprism = 0;
		int npyramid = 0;
		int ntet = 0;
		int maxlevel = 0;

		for (i = 1; i <= elements.Size(); i++)
		{
		HPRefElement el = elements.Get(i);
		maxlevel = max2(el.level, maxlevel);
		switch (Get_HPRef_Struct(el.type).geom)
		{
		  case HPREF_ELEMENT_TYPE.HP_SEGM:

		  {
			  nsegm++;
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_TRIG:
		  {
			  ntrig++;
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_QUAD:
		  {
			  nquad++;
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_TET:
		  {
			  ntet++;
			  break;
		  }

		  case HPREF_ELEMENT_TYPE.HP_PRISM:
		  {
			  nprism++;
			  break;
		  }

		  case HPREF_ELEMENT_TYPE.HP_PYRAMID:
		  {
			  npyramid++;
			  break;
		  }

		  case HPREF_ELEMENT_TYPE.HP_HEX:
		  {
			  nhex++;
			  break;
		  }

		  default:
		  {
			  cerr << "statistics error, unknown element type" << "\n";
		  }
		break;
		}
		}

		Console.Write("level = ");
		Console.Write(maxlevel);
		Console.Write("\n");
		Console.Write("nsegm = ");
		Console.Write(nsegm);
		Console.Write("\n");
		Console.Write("ntrig = ");
		Console.Write(ntrig);
		Console.Write(", nquad = ");
		Console.Write(nquad);
		Console.Write("\n");
		Console.Write("ntet = ");
		Console.Write(ntet);
		Console.Write(", npyr = ");
		Console.Write(npyramid);
		Console.Write(", nprism = ");
		Console.Write(nprism);
		Console.Write(", nhex = ");
		Console.Write(nhex);
		Console.Write("\n");

		return;

		double memcost = 0;
		double cpucost = 0;
		for (p = 1; p <= 20; p++)
		{
		memcost = (ntet + nprism + nhex) * ngsimd.GlobalMembers.pow((double)p, 6.0);
		cpucost = (ntet + nprism + nhex) * ngsimd.GlobalMembers.pow((double)p, 9.0);
		Console.Write("costs for p = ");
		Console.Write(p);
		Console.Write(": mem = ");
		Console.Write(memcost);
		Console.Write(", cpu = ");
		Console.Write(cpucost);
		Console.Write("\n");
		}

		double memcosttet = 0;
		double memcostprism = 0;
		double memcosthex = 0;
		double memcostsctet = 0;
		double memcostscprism = 0;
		double memcostschex = 0;
		double cpucosttet = 0;
		double cpucostprism = 0;
		double cpucosthex = 0;

		for (i = 1; i <= elements.Size(); i++)
		{
		HPRefElement el = elements.Get(i);
		switch (el.type)
		{
		  case HPREF_ELEMENT_TYPE.HP_TET:
		  case HPREF_ELEMENT_TYPE.HP_TET_0E_1V:
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_0V:
		  case HPREF_ELEMENT_TYPE.HP_TET_1E_1VA:
		  {
			  int p1 = maxlevel - el.level + 1;
			  (*testout) << "p1 = " << p1 << ", P1^6 = " << ngsimd.GlobalMembers.pow((double)p1, 6.0) << " (p1-3)^6 = " << ngsimd.GlobalMembers.pow((double)max2(p1 - 3, 0), 6.0) << " p1^3 = " << ngsimd.GlobalMembers.pow((double)p1, 3.0) << " (p1-3)^3 = " << ngsimd.GlobalMembers.pow((double)(p1 - 3), 3.0) << " [p1^3-(p1-3)^3]^2 = " << sqr(ngsimd.GlobalMembers.pow((double)p1, 3.0) - ngsimd.GlobalMembers.pow((double)(p1 - 3), 3.0)) << "\n";

			  p1 /= 2 + 1;
			  memcosttet += ngsimd.GlobalMembers.pow((double)p1, 6.0);
			  memcostsctet += ngsimd.GlobalMembers.pow((double)p1, 6.0) - ngsimd.GlobalMembers.pow((double)max2(p1 - 3, 1), 6.0);
			  cpucosttet += ngsimd.GlobalMembers.pow((double)p1, 9.0);
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_PRISM:
		  case HPREF_ELEMENT_TYPE.HP_PRISM_SINGEDGE:
		  {
			  int p1 = maxlevel - el.level + 1;
			  p1 /= 2 + 1;
			  memcostprism += ngsimd.GlobalMembers.pow((double)p1, 6.0);
			  memcostscprism += ngsimd.GlobalMembers.pow((double)p1, 6.0) - ngsimd.GlobalMembers.pow((double)max2(p1 - 3, 1), 6.0);
			  cpucostprism += ngsimd.GlobalMembers.pow((double)p1, 9.0);
			  break;
		  }
		  case HPREF_ELEMENT_TYPE.HP_HEX:
		  {
			  int p1 = maxlevel - el.level + 1;
			  int p2 = maxlevel;
			  p1 /= 2 + 1;
			  p2 /= 2 + 1;
			  memcosthex += ngsimd.GlobalMembers.pow((double)p1, 4.0) * ngsimd.GlobalMembers.pow((double)p2, 2.0);
			  memcostschex += ngsimd.GlobalMembers.pow((double)p1, 6.0) - ngsimd.GlobalMembers.pow((double)max2(p1 - 2, 0), 6.0);
			  cpucosthex += ngsimd.GlobalMembers.pow((double)p1, 6.0) * ngsimd.GlobalMembers.pow((double)p2, 3.0);
			  break;
		  }
		  default:
			;
		break;
		}
		}
		Console.Write("TET: hp-memcost = ");
		Console.Write(memcosttet);
		Console.Write(", scmemcost = ");
		Console.Write(memcostsctet);
		Console.Write(", cpucost = ");
		Console.Write(cpucosttet);
		Console.Write("\n");
		Console.Write("PRI: hp-memcost = ");
		Console.Write(memcostprism);
		Console.Write(", scmemcost = ");
		Console.Write(memcostscprism);
		Console.Write(", cpucost = ");
		Console.Write(cpucostprism);
		Console.Write("\n");
		Console.Write("HEX: hp-memcost = ");
		Console.Write(memcosthex);
		Console.Write(", scmemcost = ");
		Console.Write(memcostschex);
		Console.Write(", cpucost = ");
		Console.Write(cpucosthex);
		Console.Write("\n");
	#endif
	  }



	  public static void ReorderPoints(Mesh mesh, Array<HPRefElement> hpelements)
	  {
		Array<int, 1> map = new Array<int, 1>(mesh.GetNP());

		for (int i = 1; i <= mesh.GetNP(); i++)
		{
		  map[i] = i;
		}

		int nwrong = 0;
		int nright = 0;
		for (int k = 0; k < 5; k++)
		{
			nwrong = nright = 0;
			for (int i = 0; i < hpelements.Size(); i++)
			{
				HPRefElement hpel = hpelements[i];

				if (Get_HPRef_Struct(hpel.type).geom == HPREF_ELEMENT_TYPE.HP_PRISM)
				{
					int minbot = 0;
					int mintop = 0;
					for (int j = 0; j < 3; j++)
					{
						if (map[hpel.pnums[j]] < map[hpel.pnums[minbot]])
						{
							minbot = j;
						}
						if (map[hpel.pnums[j + 3]] < map[hpel.pnums[mintop + 3]])
						{
							mintop = j;
						}
					}
					if (minbot != mintop)
					{
					  nwrong++;
					}
					else
					{
					  nright++;
					}

					if (minbot != mintop)
					{
						if (map[hpel.pnums[minbot]] < map[hpel.pnums[mintop + 3]])
						{
						  swap(map[hpel.pnums[3 + minbot]], map[hpel.pnums[3 + mintop]]);
						}
						else
						{
						  swap(map[hpel.pnums[minbot]], map[hpel.pnums[mintop]]);
						}
					}
				}
			}
			// cout << nwrong << " wrong prisms, " << nright << " right prisms" << endl;
		}

		Console.Write(nwrong);
		Console.Write(" wrong prisms, ");
		Console.Write(nright);
		Console.Write(" right prisms");
		Console.Write("\n");


		Array<MeshPoint, 1> hpts = new Array<MeshPoint, 1>(mesh.GetNP());

		for (int i = 1; i <= mesh.GetNP(); i++)
		{
		  hpts[map[i]] = new mesh.Point(i);
		}

		for (int i = 1; i <= mesh.GetNP(); i++)
		{
		  mesh.Point(i) = hpts[i];
		}

		for (int i = 0; i < hpelements.Size(); i++)
		{
			HPRefElement hpel = hpelements[i];
			for (int j = 0; j < hpel.np; j++)
			{
			  hpel.pnums[j] = map[hpel.pnums[j]];
			}
		}
	  }




	  public static ostream operator << (ostream ost, GradingBox box)
	  {
		ost << "gradbox, pmid = " << box.PMid() << ", h2 = " << box.H2() << " cutbound = " << box.flags.cutboundary << " isinner = " << box.flags.isinner << "\n";
		return ost;
	  }



	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int FindInnerBoxes_timer = NgProfiler.CreateTimer("LocalH::FindInnerBoxes");
















	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int FindInnerBoxes_timer = NgProfiler.CreateTimer("LocalH::FindInnerBoxes 2d");

	  public static ostream operator << (ostream ost, Mesh mesh)
	  {
		ost << "mesh: " << "\n";
		mesh.Save(ost);
		return ost;
	  }

	  internal static object buildsearchtree_mutex = new object();


	  public static int SolveLinearSystemLS(Vec3d col1, Vec3d col2, Vec3d rhs, ref Vec2d sol)
	  {
		double a11 = col1 * col1;
		double a12 = col1 * col2;
		double a22 = col2 * col2;

		double det = a11 * a22 - a12 * a12;

		if (det * det <= 1e-24 * a11 * a22)
		{
			sol = new Vec2d(0, 0);
			return 1;
		}

		Vec2d aTrhs = new Vec2d();
		aTrhs.X() = col1 * rhs;
		aTrhs.Y() = col2 * rhs;

		sol.X() = (a22 * aTrhs.X() - a12 * aTrhs.Y()) / det;
		sol.Y() = (-a12 * aTrhs.X() + a11 * aTrhs.Y()) / det;
		return 0;
	  }

	  public static bool ValidBarCoord(double[] lami, double eps = 1e-12)
	  {
		return (lami[0] <= 1.0 + eps && lami[0] >= 0.0 - eps && lami[1] <= 1.0 + eps && lami[1] >= 0.0 - eps && lami[2] <= 1.0 + eps && lami[2] >= 0.0 - eps);
	  }

	  public static string Mesh.defaultmat = "default";

	  public static string Mesh.cd2_default_name = "default";
	  public static string Mesh.default_bc = "default";

	  public static string Mesh.cd3_default_name = "default";
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern const char * tetrules[];
	  // extern const char * tetrules2[];
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern const char * prismrules2[];
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern const char * pyramidrules[];
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern const char * pyramidrules2[];
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern const char * hexrules[];

	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  internal static int Optimize2d_timer = NgProfiler.CreateTimer("optimize2d");

	  public static DLL_HEADER void Optimize2d(Mesh mesh, MeshingParameters mp)
	  {
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//	static int timer = NgProfiler::CreateTimer("optimize2d");
		NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(Optimize2d_timer);

		mesh.CalcSurfacesOfNode();

		char[] optstr = mp.optimize2d;
		int optsteps = mp.optsteps2d;

		for (int i = 1; i <= optsteps; i++)
		{
		  for (uint j = 1; j <= optstr.Length; j++)
		  {
		  if (multithread.terminate)
		  {
			  break;
		  }
		  switch (optstr[j - 1])
		  {
			case 's':
			{ // topological swap
			MeshOptimize2d meshopt = new MeshOptimize2d();
					meshopt.SetMetricWeight(mp.elsizeweight);
			meshopt.EdgeSwapping(mesh, 0);
			break;
			}
			case 'S':
			{ // metric swap
			MeshOptimize2d meshopt = new MeshOptimize2d();
					meshopt.SetMetricWeight(mp.elsizeweight);
			meshopt.EdgeSwapping(mesh, 1);
			break;
			}
			case 'm':
			{
			MeshOptimize2d meshopt = new MeshOptimize2d();
					meshopt.SetMetricWeight(mp.elsizeweight);
			meshopt.ImproveMesh(mesh, mp);
			break;
			}
			case 'c':
			{
			MeshOptimize2d meshopt = new MeshOptimize2d();
					meshopt.SetMetricWeight(mp.elsizeweight);
			meshopt.CombineImprove(mesh);
			break;
			}
			default:
			  cerr << "Optimization code " << optstr[j - 1] << " not defined" << "\n";
		  break;
		  }
		  }
		}
	  }
	  internal static void glrender(int wait)
	  {
		//  cout << "plot adfront" << endl;

		if (multithread.drawing)
		{
		//      vssurfacemeshing.Render();
		Render();

		if (wait != 0 || multithread.testmode)
		{
			multithread.pause = 1;
		}
//C++ TO C# CONVERTER TODO TASK: Variables cannot be declared in if/while/switch conditions in C#:
		while (multithread.pause)
		{
			;
		}
		}
	  }

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern STLGeometry * stlgeometry;
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern Mesh * mesh;
	  public static VisualSceneSurfaceMeshing vssurfacemeshing = new VisualSceneSurfaceMeshing();
	  public static void glrender(int wait)
	  {
		  ;
	  }

	public static double minother;
	public static double minwithoutother;





	  internal static double TriangleQualityInst(Point3d p1, Point3d p2, Point3d p3)
	  {
		// quality 0 (worst) .. 1 (optimal)

		Vec3d v1 = new Vec3d();
		Vec3d v2 = new Vec3d();
		Vec3d v3 = new Vec3d();
		double s1;
		double s2;
		double s3;
		double an1;
		double an2;
		double an3;

		v1 = p2 - p1;
		v2 = p3 - p1;
		v3 = p3 - p2;

		an1 = Angle(v1, v2);
		v1 *= -1;
		an2 = Angle(v1, v3);
		an3 = Angle(v2, v3);

		s1 = ngsimd.GlobalMembers.sin(an1 / 2);
		s2 = ngsimd.GlobalMembers.sin(an2 / 2);
		s3 = ngsimd.GlobalMembers.sin(an3 / 2);

		return 8 * s1 * s2 * s3;
	  }


	  internal static double TetElementQuality(Point3d p1, Point3d p2, Point3d p3, Point3d p4)
	  {
		double vol;
		double l;
		double l4;
		double l5;
		double l6;


		Vec3d v1 = p2 - p1;
		Vec3d v2 = p3 - p1;
		Vec3d v3 = p4 - p1;

		vol = ngsimd.GlobalMembers.fabs((Cross(new netgen.Vec3d(v1), new netgen.Vec3d(v2)) * v3)) / 6;
		l4 = Dist(p2, p3);
		l5 = Dist(p2, p4);
		l6 = Dist(p3, p4);

		l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

		if (vol <= 1e-8 * l * l * l)
		{
			return 1e-10;
		}

		return vol / (l * l * l) * 1832.82; // 6^4 * sqrt(2)
	  }


	#if OLD
	  public static void Save2DMesh(Mesh mesh2d, Array<SplineSegment > splines, ostream outfile)

	  {
		int i;
		int j;
		outfile.precision(6);

		outfile << "areamesh2" << "\n";


		outfile << "\n";
		outfile << mesh2d.GetNSeg() << "\n";
		for (i = 1; i <= mesh2d.GetNSeg(); i++)
		{
		  outfile << mesh2d.LineSegment(i).si << "        " << mesh2d.LineSegment(i)[0] << " " << mesh2d.LineSegment(i)[1] << "  " << "\n";
		}


		outfile << mesh2d.GetNSE() << "\n";
		for (i = 1; i <= mesh2d.GetNSE(); i++)
		{
		outfile << mesh2d.SurfaceElement(i).GetIndex() << "         ";
		outfile << mesh2d.SurfaceElement(i).GetNP() << " ";
		for (j = 1; j <= mesh2d.SurfaceElement(i).GetNP(); j++)
		{
		  outfile << mesh2d.SurfaceElement(i).PNum(j) << " ";
		}
		outfile << "\n";
		}

		outfile << mesh2d.GetNP() << "\n";
		for (i = 1; i <= mesh2d.GetNP(); i++)
		{
		  outfile << new mesh2d.Point(i).X() << " " << new mesh2d.Point(i).Y() << "\n";
		}

		if (splines != null)
		{
		outfile << splines.Size() << "\n";
		for (i = 1; i <= splines.Size(); i++)
		{
		  splines.Get(i).PrintCoeff(outfile);
		}
		}
		else
		{
		  outfile << "0" << "\n";
		}
	  }
	#endif








	  public static void SaveVolumeMesh(Mesh mesh, NetgenGeometry geometry, ref string filename)
	  {
		int i;

		ofstream outfile = new ofstream(filename);
		outfile << "volumemesh" << "\n";

		outfile << mesh.GetNSE() << "\n";
		for (i = 1; i <= mesh.GetNSE(); i++)
		{
		if ((mesh.SurfaceElement(i).GetIndex()) != 0)
		{
		  outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex()).SurfNr() << "\t";
		}
		else
		{
		  outfile << "0" << "\t";
		}
		outfile << mesh.SurfaceElement(i)[0] << " " << mesh.SurfaceElement(i)[1] << " " << mesh.SurfaceElement(i)[2] << "\n";
		}
		outfile << mesh.GetNE() << "\n";
		for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
		{
		  outfile << mesh[ei].GetIndex() << "\t" << mesh[ei][0] << " " << mesh[ei][1] << " " << mesh[ei][2] << " " << mesh[ei][3] << "\n";
		}

		outfile << mesh.GetNP() << "\n";
		for (i = 1; i <= mesh.GetNP(); i++)
		{
		  outfile << new mesh.Point(i)(0) << " " << new mesh.Point(i)(1) << " " << new mesh.Point(i)(2) << "\n";
		}

	#if SOLIDGEOM
		outfile << geometry.GetNSurf() << "\n";
		for (i = 1; i <= geometry.GetNSurf(); i++)
		{
		  geometry.GetSurface(i).Print(outfile);
		}
	#endif
	  }


//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern DLL_HEADER size_t timestamp;
	  public static uint GetTimeStamp()
	  {
		return timestamp;
	  }

	  public static uint NextTimeStamp()
	  {
		timestamp++;
		return timestamp;
	  }

	  public static ostream operator << (ostream ost, PointGeomInfo gi)
	  {
		return (ost << gi.trignum << " " << gi.u << " " << gi.v);
	  }

	  public static istream operator >> (istream ist, PointGeomInfo gi)
	  {
		return (ist >> gi.trignum >> gi.u >> gi.v);
	  }

	  public static ostream operator << (ostream ost, EdgePointGeomInfo gi)
	  {
		ost << "epgi: edgnr=" << gi.edgenr << ", dist=" << gi.dist;
		return ost;
	  }

	  public static istream operator >> (istream ist, ref PointIndex pi)
	  {
		int i;
		ist >> i;
		pi = new PointIndex(i);
		return ist;
	  }

	  public static ostream operator << (ostream ost, PointIndex pi)
	  {
		return (ost << (int)pi);
	  }

//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <int N>
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
	//  class PointIndices;

	  public static istream operator >> (istream ist, ref ElementIndex pi)
	  {
		int i;
		ist >> i;
		pi = i;
		return ist;
	  }

	  public static ostream operator << (ostream ost, ElementIndex pi)
	  {
		return (ost << (int)pi);
	  }

	  public static istream operator >> (istream ist, ref SurfaceElementIndex pi)
	  {
		int i;
		ist >> i;
		pi = i;
		return ist;
	  }

	  public static ostream operator << (ostream ost, SurfaceElementIndex pi)
	  {
		return (ost << (int)pi);
	  }

	  public static istream operator >> (istream ist, ref SegmentIndex pi)
	  {
		int i;
		ist >> i;
		pi = i;
		return ist;
	  }

	  public static ostream operator << (ostream ost, SegmentIndex pi)
	  {
		return (ost << (int)pi);
	  }

	  public static ostream operator << (ostream s, MeshPoint pt)
	  {
		return (s << Point < 3> (pt));
	  }

	  public static ostream operator << (ostream s, Element2d el)
	  {
		s << "np = " << el.GetNP();
		for (int j = 1; j <= el.GetNP(); j++)
		{
		  s << " " << el.PNum(j);
		}
		return s;
	  }

	  public static ostream operator << (ostream s, Element el)
	  {
		s << "np = " << el.GetNP();
		for (int j = 0; j < el.GetNP(); j++)
		{
		  s << " " << (int)el[j];
		}
		return s;
	  }

	  public static ostream operator << (ostream s, Segment seg)
	  {
		s << seg[0] << "(gi=" << seg.geominfo[0].trignum << ") - " << seg[1] << "(gi=" << seg.geominfo[1].trignum << ")" << " domin = " << seg.domin << ", domout = " << seg.domout << " si = " << seg.si << ", edgenr = " << seg.edgenr;
		return s;
	  }

	  public static ostream operator << (ostream s, Element0d el)
	  {
		s << el.pnum << ", index = " << el.index;
		return s;
	  }

	  public static ostream operator << (ostream s, FaceDescriptor fd)
	  {
		s << "surfnr = " << fd.SurfNr() << ", domin = " << fd.DomainIn() << ", domout = " << fd.DomainOut() << ", tlosurf = " << fd.TLOSurface() << ", bcprop = " << fd.BCProperty() << ", bcname = " << fd.GetBCName() << ", domin_sing = " << fd.DomainInSingular() << ", domout_sing = " << fd.DomainOutSingular() << ", colour = " << fd.SurfColour();
		return s;
	  }

	  public static ostream operator << (ostream ost, MeshingParameters mp)
	  {
		mp.Print(ost);
		return ost;
	  }



	  internal int[][] gftetfacesa =
	  {
		  new int[] {1, 2, 3},
		  new int[] {2, 0, 3},
		  new int[] {0, 1, 3},
		  new int[] {1, 0, 2}
	  };







	  public static Array<IntegrationPointData> ipdtrig = new Array<IntegrationPointData>();
	  public static Array<IntegrationPointData> ipdquad = new Array<IntegrationPointData>();



	  internal int[][] qip_table =
	  {
		  new int[] {0, 1, 0, 3},
		  new int[] {0, 1, 1, 2},
		  new int[] {3, 2, 0, 3},
		  new int[] {3, 2, 1, 2}
	  };





	  public static Array< IntegrationPointData> ipdtet = new Array< IntegrationPointData>();
	  public static Array< IntegrationPointData> ipdtet10 = new Array< IntegrationPointData>();

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: template void Element2d::GetShapeNew(const Point<2,double> & p, TFlatVector<double> shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element2d::GetShapeNew(Point<2,double> p, TFlatVector<double> shape);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: template void Element2d::GetShapeNew(const Point<2,SIMD<double>> & p, TFlatVector<SIMD<double>> shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element2d::GetShapeNew(Point<2,SIMD<double>> p, TFlatVector<SIMD<double>> shape);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void Element2d::GetDShapeNew<double> (const Point<2> &, MatrixFixWidth<2> &) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element2d::GetDShapeNew (Point<2> UnnamedParameter, MatrixFixWidth<2> UnnamedParameter2);
//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void Element2d::GetDShapeNew<SIMD<double>> (const Point<2,SIMD<double>> &, MatrixFixWidth<2,SIMD<double>> &) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element2d::GetDShapeNew (Point<2,SIMD<double>> UnnamedParameter, MatrixFixWidth<2,SIMD<double>> UnnamedParameter2);


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: template void Element::GetShapeNew(const Point<3,double> & p, TFlatVector<double> shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element::GetShapeNew(Point<3,double> p, TFlatVector<double> shape);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: template void Element::GetShapeNew(const Point<3,SIMD<double>> & p, TFlatVector<SIMD<double>> shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element::GetShapeNew(Point<3,SIMD<double>> p, TFlatVector<SIMD<double>> shape);

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void Element::GetDShapeNew<double> (const Point<3> &, MatrixFixWidth<3> &) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element::GetDShapeNew (Point<3> UnnamedParameter, MatrixFixWidth<3> UnnamedParameter2);
//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: template void Element::GetDShapeNew<SIMD<double>> (const Point<3,SIMD<double>> &, MatrixFixWidth<3,SIMD<double>> &) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//  template void Element::GetDShapeNew (Point<3,SIMD<double>> UnnamedParameter, MatrixFixWidth<3,SIMD<double>> UnnamedParameter2);

  //the dots for progression of program

	  public static void PrintDot(char ch = '.')
	  {
		// if (printdots)
		if (printmessage_importance >= 4)
		{
			string st = new string(new char[2]);
			st = StringFunctions.ChangeCharacter(st, 0, ch);
			st = st.Substring(0, 1);
			Ng_PrintDest(st);
		}
	  }


	  //Message Pipeline:

	  //importance: importance of message: 1=very important, 3=middle, 5=low, 7=unimportant
	  public static void PrintMessage(int importance, MyStr s1, MyStr s2 = new MyStr())
	  {
		if (importance <= printmessage_importance)
		{
			Ng_PrintDest(new MyStr(" ") + s1.functorMethod + s2.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintMessage(int importance, MyStr s1, MyStr s2, MyStr s3, MyStr s4 = new MyStr())
	  {
		if (importance <= printmessage_importance)
		{
			Ng_PrintDest(new MyStr(" ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintMessage(int importance, MyStr s1, MyStr s2, MyStr s3, MyStr s4, MyStr s5, MyStr s6 = new MyStr(), MyStr s7 = new MyStr(), MyStr s8 = new MyStr())
	  {
		if (importance <= printmessage_importance)
		{
			Ng_PrintDest(new MyStr(" ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

	  // CR without line-feed
	  public static void PrintMessageCR(int importance, MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (importance <= printmessage_importance)
		{
			Ng_PrintDest(new MyStr(" ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\r"));
		}
	  }

	  public static void PrintFnStart(MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (printfnstart)
		{
		  Ng_PrintDest(new MyStr(" Start Function: ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintWarning(MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (printwarnings)
		{
		  Ng_PrintDest(new MyStr(" WARNING: ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintError(MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (printerrors)
		{
		  Ng_PrintDest(new MyStr(" ERROR: ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintFileError(MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (printerrors)
		{
		  Ng_PrintDest(new MyStr(" FILE ERROR: ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintSysError(MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (printerrors)
		{
		  Ng_PrintDest(new MyStr(" SYSTEM ERROR: ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

	  public static void PrintUserError(MyStr s1, MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		Ng_PrintDest(new MyStr(" USER ERROR: ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
	  }

	  public static void PrintTime(MyStr s1 = "", MyStr s2 = "", MyStr s3 = "", MyStr s4 = "", MyStr s5 = "", MyStr s6 = "", MyStr s7 = "", MyStr s8 = "")
	  {
		if (printmessage_importance >= 3)
		{
		  Ng_PrintDest(new MyStr(" Time = ") + s1.functorMethod + s2.functorMethod + s3.functorMethod + s4.functorMethod + s5.functorMethod + s6.functorMethod + s7.functorMethod + s8.functorMethod + new MyStr("\n"));
		}
	  }

  /*
  void SetStatMsgF(const MyStr& s)
  {
    PrintFnStart(s);
    SetStatMsg(s);
  }
  */


	  public static void SetStatMsg(MyStr s)
	  {
		msgstatus.CopyFrom(s.functorMethod);
		multithread.task = msgstatus.c_str();
	  }

	  public static void PushStatus(MyStr s)
	  {
		msgstatus_stack.Append(new MyStr(s.functorMethod));
		SetStatMsg(s.functorMethod);
		threadpercent_stack.Append(0);
	  }

	  public static void PushStatusF(MyStr s)
	  {
		msgstatus_stack.Append(new MyStr(s.functorMethod));
		SetStatMsg(s.functorMethod);
		threadpercent_stack.Append(0);
		PrintFnStart(s.functorMethod);
	  }

	  public static void PopStatus()
	  {
		if (msgstatus_stack.Size())
		{
			if (msgstatus_stack.Size() > 1)
			{
		  // SetStatMsg (*msgstatus_stack.Last());
		  SetStatMsg(msgstatus_stack[msgstatus_stack.Size() - 2]);
			}
			else
			{
		  SetStatMsg("");
			}
			msgstatus_stack.Last() = null;
			msgstatus_stack.DeleteLast();
			threadpercent_stack.DeleteLast();
			if (threadpercent_stack.Size() > 0)
			{
		  multithread.percent = threadpercent_stack.Last();
			}
			else
			{
		  multithread.percent = 100.0;
			}
		}
		else
		{
			PrintSysError("PopStatus failed");
		}
	  }

	  public static void SetThreadPercent(double percent)
	  {
		multithread.percent = percent;
		if (threadpercent_stack.Size() > 0)
		{
		  threadpercent_stack.Last() = percent;
		}
	  }

	  public static void GetStatus(ref MyStr s, ref double percentage)
	  {
		if (threadpercent_stack.Size() > 0)
		{
		  percentage = threadpercent_stack.Last();
		}
		else
		{
		  percentage = multithread.percent;
		}

		if (msgstatus_stack.Size())
		{
		  s = *msgstatus_stack.Last();
		}
		else
		{
		  s = "idle";
		}
	  }
	public static int printmessage_importance = 5;
	public static int printwarnings = 1;
	public static int printerrors = 1;
	public static int printdots = 1;
	public static int printfnstart = 0;

	// extern void Ng_PrintDest(const MyStr& s);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void Ng_PrintDest(string s);


//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	static Array<MyStr*> msgstatus_stack(0);
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	static Array<double> threadpercent_stack(0);
	internal static MyStr msgstatus = "";




	public static void ResetStatus()
	{
	  SetStatMsg("idle");

	  for (int i = 0; i < msgstatus_stack.Size(); i++)
	  {
		msgstatus_stack[i] = null;
	  }
	  msgstatus_stack.SetSize(0);
	  threadpercent_stack.SetSize(0);

	  // multithread.task = "";
	  multithread.percent = 100.0;
	}

//C++ TO C# CONVERTER TODO TASK: C++ template specialization was removed by C++ to C# Converter:
//ORIGINAL LINE: inline int MyGetMPIType<PointIndex> ()
	  public static int MyGetMPIType ()
	  {
		  return MPI_INT;
	  }


	public static void LoadMatrixLine(istream ist, DenseMatrix m, int line)
	{
	  char ch;
	  int pnum;
	  float f;

	  ist >> ch;
	  while (ch != '}')
	  {
		  ist.putback(ch);
		  ist >> f;
		  ist >> ch;
		  ist >> pnum;

		  if (ch == 'x' || ch == 'X')
		  {
		m.Elem(line, 2 * pnum - 1) = f;
		  }
		  if (ch == 'y' || ch == 'Y')
		  {
		m.Elem(line, 2 * pnum) = f;
		  }

		  ist >> ch;
		  if (ch == DefineConstants.COMMASIGN)
		  {
		ist >> ch;
		  }
	  }
	}




//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//extern const char * triarules[];
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//extern const char * quadrules[];

//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//extern const char * tetrules[];

	public static void LoadVMatrixLine(istream ist, DenseMatrix m, int line)
	{
	  char ch;
	  int pnum;
	  float f;

	  ist >> ch;
	  while (ch != '}')
	  {
		  ist.putback(ch);
		  ist >> f;
		  ist >> ch;
		  ist >> pnum;

		  if (ch == 'x' || ch == 'X')
		  {
		m.Elem(line, 3 * pnum - 2) = f;
		  }
		  if (ch == 'y' || ch == 'Y')
		  {
		m.Elem(line, 3 * pnum - 1) = f;
		  }
		  if (ch == 'z' || ch == 'Z')
		  {
		m.Elem(line, 3 * pnum) = f;
		  }

		  if (ch == 'p' || ch == 'P')
		  {
		  m.Elem(line, 3 * pnum - 2) = f;
		  m.Elem(line+1, 3 * pnum - 1) = f;
		  m.Elem(line+2, 3 * pnum) = f;
		  }

		  ist >> ch;
		  if (ch == DefineConstants.COMMASIGN)
		  {
		ist >> ch;
		  }
	  }
	}
	public static string[] prismrules2 = {"tolfak 0.5\n", \ "\n", \ "rule \"prism on quad\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3, 6);\n", \ "(1, 5, 6, 4);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 5);\n", \ "(1, 3, 4, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "{ -0.1 P1, 0.3 P2, 0.3 P3, -0.1 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "{ 0 P1, 0.25 P2, 0.25 P3, 0 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5 6 7;\n", \ "\n", \ "freeset\n", \ "2 3 4 5 6 8;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on quad, one trig\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(1, 5, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3, 6);\n", \ "(1, 5, 6, 4);\n", \ "(4, 6, 3);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 5);\n", \ "(1, 3, 4, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "{ -0.1 P1, 0.3 P2, 0.3 P3, -0.1 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "{ 0 P1, 0.25 P2, 0.25 P3, 0 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5 6 7;\n", \ "\n", \ "freeset\n", \ "2 3 4 5 6 8;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 2 quad\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 5, 6, 4);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5 6 7;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 2 quad, one trig\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(1, 5, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 5, 6, 4);\n", \ "(4, 6, 3);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5 6 7;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 2 quada\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(5, 1, 4, 6) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3, 6);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ -0.1 P1, 0.3 P2, 0.3 P3, -0.1 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0 P1, 0.25 P2, 0.25 P3, 0 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5 6 7;\n", \ "\n", \ "freeset\n", \ "1 3 4 6;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"fill prism\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(5, 1, 4, 6) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 3 quad, one trig\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(5, 1, 4, 6) del;\n", \ "(1, 5, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(4, 6, 3);\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"flat prism\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(0.5, 0.866, 0);\n", \ "(0, 0, -1);\n", \ "(1, 0, -1);\n", \ "(0.5, 0.866, -1);\n", \ "\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(5, 4, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 2, 4);\n", \ "(4, 2, 5);\n", \ "(2, 3, 5);\n", \ "(5, 3, 6);\n", \ "(3, 1, 6);\n", \ "(6, 1, 4);\n", \ "\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 5, 4, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "endrule\n", \ "\n", \ 0};
	public static string[] prismrules2 = {"tolfak 0.5\n", \ "\n", \ "rule \"prism on quad\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3, 6);\n", \ "(1, 5, 6, 4);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 5);\n", \ "(1, 3, 4, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "{ -0.1 P1, 0.3 P2, 0.3 P3, -0.1 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "{ 0 P1, 0.25 P2, 0.25 P3, 0 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on quad, one trig\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(1, 5, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3, 6);\n", \ "(1, 5, 6, 4);\n", \ "(4, 6, 3);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 5);\n", \ "(1, 3, 4, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "{ -0.1 P1, 0.3 P2, 0.3 P3, -0.1 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "{ 0 P1, 0.25 P2, 0.25 P3, 0 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 2 quad\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 5, 6, 4);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5 6 7;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 2 quad, one trig\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(1, 5, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 5, 6, 4);\n", \ "(4, 6, 3);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.3 P1, -0.1 P2, -0.1 P3, 0.3 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0.25 P1, 0 P2, 0 P3, 0.25 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5 6 7;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 2 quada\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(5, 1, 4, 6) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3, 6);\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ -0.1 P1, 0.3 P2, 0.3 P3, -0.1 P4, 0.3 P5, 0.3 P6 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 0 P1, 0.25 P2, 0.25 P3, 0 P4, 0.25 P5, 0.25 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5 6 7;\n", \ "\n", \ "freeset\n", \ "1 3 4 6;\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"fill prism\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(5, 1, 4, 6) del;\n", \ "(1, 5, 2) del;\n", \ "(4, 3, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"prism on 3 quad, one trig\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0, -0.86);\n", \ "(0.5, 1, -0.86);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 5, 6, 3) del;\n", \ "(5, 1, 4, 6) del;\n", \ "(1, 5, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(4, 6, 3);\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 5, 2, 4, 6, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "freeset\n", \ "1 2 4 5;\n", \ "\n", \ "freeset\n", \ "2 3 4 6;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"flat prism\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(0.5, 0.866, 0);\n", \ "(0, 0, -1);\n", \ "(1, 0, -1);\n", \ "(0.5, 0.866, -1);\n", \ "\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(5, 4, 6) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 2, 4);\n", \ "(4, 2, 5);\n", \ "(2, 3, 5);\n", \ "(5, 3, 6);\n", \ "(3, 1, 6);\n", \ "(6, 1, 4);\n", \ "\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 5, 4, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "endrule\n", \ "\n", \ 0};
	public static string[] pyramidrules2 = {"tolfak 0.5\n", \ "\n", \ "rule \"Pyramid on quad\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.5, -0.5) \n", \ "	{ 0.25 X1, 0.25 X2, 0.25 X3, 0.25 X4 } 	\n", \ "	{ 0.25 Y1, 0.25 Y2, 0.25 Y3, 0.25 Y4 } { };\n", \ "\n", \ "newfaces\n", \ "(1, 2, 5);\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.4 P5, -0.1 P1, -0.1 P2, -0.1 P3, -0.1 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"small Pyramid on quad\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.5, -0.1 )\n", \ "	{ 0.25 X1, 0.25 X2, 0.25 X3, 0.25 X4 } \n", \ "	{ 0.25 Y1, 0.25 Y2, 0.25 Y3, 0.25 Y4 } { };\n", \ "\n", \ "newfaces\n", \ "(1, 2, 5);\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.4 P5, -0.1 P1, -0.1 P2, -0.1 P3, -0.1 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"connect pyramid\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 2, 5);\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"pyramid with one trig\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 1, 5) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 0.34 P2, 0.34 P3, 0.34 P5, -0.02 P1 };\n", \ "{ 0.34 P3, 0.34 P4, 0.34 P5, -0.02 P1 };\n", \ "{ 0.34 P1, 0.34 P4, 0.34 P5, -0.02 P2 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 0.333 P2, 0.333 P3, 0.334 P5, 0 P1 };\n", \ "{ 0.333 P3, 0.333 P4, 0.334 P5, 0 P1 };\n", \ "{ 0.333 P1, 0.333 P4, 0.334 P5, 0 P2 };\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 5);\n", \ "(1, 3, 4, 5);\n", \ "\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "freeset\n", \ "1 3 4 5;\n", \ "freeset\n", \ "2 3 5 6;\n", \ "freeset\n", \ "3 4 5 7;\n", \ "freeset \n", \ "1 4 5 8;\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"pyramid with two trig\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 1, 5) del;\n", \ "(3, 2, 5) del;\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"pyramid with two trig, left\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 1, 5) del;\n", \ "(1, 4, 5) del;\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(3, 4, 5);\n", \ "(2, 3, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ 0};
	public static string[] pyramidrules = {"tolfak 0.5\n", \ "\n", \ "rule \"Pyramid on quad\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.5, -0.5) \n", \ "	{ 0.25 X1, 0.25 X2, 0.25 X3, 0.25 X4 } 	\n", \ "	{ 0.25 Y1, 0.25 Y2, 0.25 Y3, 0.25 Y4 } { };\n", \ "\n", \ "newfaces\n", \ "(1, 2, 5);\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.4 P5, -0.1 P1, -0.1 P2, -0.1 P3, -0.1 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"small Pyramid on quad\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.5, -0.1 )\n", \ "	{ 0.25 X1, 0.25 X2, 0.25 X3, 0.25 X4 } \n", \ "	{ 0.25 Y1, 0.25 Y2, 0.25 Y3, 0.25 Y4 } { };\n", \ "\n", \ "newfaces\n", \ "(1, 2, 5);\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.4 P5, -0.1 P1, -0.1 P2, -0.1 P3, -0.1 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"connect pyramid\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 2, 5);\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"pyramid with one trig\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 1, 5) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(2, 3, 5);\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 0.34 P2, 0.34 P3, 0.34 P5, -0.02 P1 };\n", \ "{ 0.34 P3, 0.34 P4, 0.34 P5, -0.02 P1 };\n", \ "{ 0.34 P1, 0.34 P4, 0.34 P5, -0.02 P3 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 0.333 P2, 0.333 P3, 0.334 P5, 0 P1 };\n", \ "{ 0.333 P3, 0.333 P4, 0.334 P5, 0 P1 };\n", \ "{ 0.333 P1, 0.333 P4, 0.334 P5, 0 P3 };\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 5);\n", \ "(1, 3, 4, 5);\n", \ "\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "freeset\n", \ "1 3 4 5;\n", \ "freeset\n", \ "2 3 5 6;\n", \ "freeset\n", \ "3 4 5 7;\n", \ "freeset \n", \ "1 4 5 8;\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"pyramid with two trig\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(1, 1, 0);\n", \ "(0, 1, 0);\n", \ "(0.5, 0.5, -0.5);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3, 4) del;\n", \ "(2, 1, 5) del;\n", \ "(3, 2, 5) del;\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(3, 4, 5);\n", \ "(4, 1, 5);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ 0};
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern bool netgen_executable_started;
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//  extern NetgenGeometry *ng_geometry;
	#if PARALLEL
	  /** we need allreduce in python-wrapped communicators **/
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  public static T MyMPI_AllReduceNG<T>(T d, MPI_Op op, MPI_Comm comm)
	  {
		T global_d = new default(T);
		MPI_Allreduce(d, global_d, 1, MyGetMPIType<T>(), op, comm);
		return global_d;
	  }
	#else
	  // enum { MPI_SUM = 0, MPI_MIN = 1, MPI_MAX = 2 };
	  // typedef int MPI_Op;
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
	  public static T MyMPI_AllReduceNG<T>(T d, MPI_Op op, MPI_Comm comm)
	  {
		  return d;
	  }
	#endif
	public static string[] quadrules = {"rule \"Free Quad (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(1, 1) { 1 X2 } { };\n", \ "(0, 1) { } { };\n", \ "\n", \ "newlines\n", \ "(3, 2);\n", \ "(4, 3);\n", \ "(1, 4);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 1.5) { 1.5 X2 } { };\n", \ "(-0.5, 1.5) { -0.5 X2 } { };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Free Quad (5)\"\n", \ "\n", \ "quality 5\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(1, 1) { 1 X2 } { };\n", \ "(0, 1) { } { };\n", \ "\n", \ "newlines\n", \ "(3, 2);\n", \ "(4, 3);\n", \ "(1, 4);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 1.5) { 1.5 X2 } { };\n", \ "(-0.5, 1.5) { -0.5 X2 } { };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X2 } { };\n", \ "(0, 1) { } { };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Quad Right (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "\n", \ "newpoints\n", \ "(0, 1) { } { 1 y3 };\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "(4, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(-0.5, 1.5) { } { 1.5 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Quad P Right (2)\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0, 1) { -1 X2, 1 X3 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "(4, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.2, 0.5) { 0.7 X2, 0.5 X3 } { 0.5 Y3 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(-0.5, 1.5) { -2 X2, 1.5 X3 } { 1.5 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { -1 X2, 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "\n", \ "orientations\n", \ "(1, 2, 3);\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"Quad P Right (150)\"\n", \ "\n", \ "quality 150\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0, 1) { 1 X2, -1 X3 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(1, 4)\n;", \ "(4, 3)\n;", \ "(3, 2)\n;", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.2, 0.5) { 0.7 X2, 0.5 X3 } { 0.5 Y3 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(-0.5, 1.5) { -2 X2, 1.5 X3 } { 1.5 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X2, -1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "orientations\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"Quad Right PL (2)\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "(4, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0.5, 1.2) { -0.1 X2, 0.6 X3, 0.6 X4 } { -0.1 Y2, 0.6 Y3, 0.6 Y4 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "(-0.2, 0.5) { -0.1 X2, -0.1 X3, 0.6 X4 } { -0.1 Y2, -0.1 Y3, 0.6 Y4 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0.5, 1) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "(0, 0.5) { 0.5 X4 } { 0.5 Y4 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "orientations\n", \ "(1, 2, 3);\n", \ "(1, 3, 4);\n", \ "(1, 2, 4);\n", \ "(4, 2, 3);\n", \ "\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left Quad (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(3, 4);\n", \ "(4, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 1.5) { 1.5 X2, 1.5 X3 } { 1.5 Y3 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left P Quad (2)\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 4);\n", \ "(4, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 1.5) { 1.5 X2, 1.5 X3 } { 1.5 Y3 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "(-0.2, 0.6) { -0.2 X2, 0.6 X3 } { 0.6 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 0.5) { 0.5 X3 } { 0.5 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left P Quad (150)\"\n", \ "\n", \ "quality 150\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(1, 1) { 1 X2, -1 X3 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 4);\n", \ "(4, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 1.5) { 1.5 X2, 1.5 X3 } { 1.5 Y3 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "(-0.2, 0.6) { -0.2 X2, 0.6 X3 } { 0.6 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X2, -1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 0.5) { 0.5 X3 } { 0.5 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left Quad RP (2)\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0, 1);\n", \ "(1, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(3, 4);\n", \ "(4, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.2, 0.5) { 0.6 X2, 0.6 X4, -0.1 X3 } { 0.6 Y2, 0.6 Y4, -0.1 Y3 };\n", \ "(1, 1) { 1 X4 } { 1 Y4 };\n", \ "(0.5, 1.2) { -0.1 X2, 0.6 X3, 0.6 X4 } { -0.1 Y2, 0.6 Y3, 0.6 Y4 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.5) { 0.5 X2, 0.5 X4 } { 0.5 Y2, 0.5 Y4 };\n", \ "(1, 1) { 1 X4 } { 1 Y4 };\n", \ "(0.5, 1) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n", \ "(0, 1) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4, 3);\n", \ "\n", \ "orientations\n", \ "(1, 2, 4);\n", \ "(1, 4, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Two left (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 4) del;\n", \ "(4, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.5) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y3, -0.25 Y4 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Two Right (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "(3, 4) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "(-0.5, 0.5) { -0.25 X2, -0.25 X3, 0.75 X4 } { -0.25 Y3, 0.75 Y4 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "(0, 0.5) { 0.5 X4 } { 0.5 Y4 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Right 120 (1)\"\n", \ "\n", \ "quality 1000\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 1 X3, -1 X2 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "(4, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(1, 1.732) { -2 X2, 2 X3 } { 2 Y3 };\n", \ "(0, 1.732) { -3 X2, 2 X3 } { 2 Y3 };\n", \ "(-0.5, 0.866) { -2 X2, 1 X3 } {1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4);\n", \ "(2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left 120 (1)\"\n", \ "\n", \ "quality 1000\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(-0.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 1 X3, 1 X2 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(3, 4);\n", \ "(4, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.866) { 2 X2, 1 X3 } { 1 Y3 };\n", \ "(1, 1.732) { 2 X2, 2 X3 } { 2 Y3 };\n", \ "(0, 1.732) { -1 X2, 2 X3 } { 2 Y3 };\n", \ "(-0.5, 0.866) { 1 X3 } {1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4);\n", \ "(2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left Right (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "(4, 1) del;\n", \ "\n", \ "\n", \ "newlines\n", \ "(4, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0.5, 1.5) { -0.25 X2, 0.75 X3, 0.75 X4 } { 0.75 Y3, 0.75 Y4 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0.5, 1) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Fill Quad\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "(3, 4) del;\n", \ "(4, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { 1 Y2 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Fill Triangle\"\n", \ "\n", \ "quality 10\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0.5, 0.86);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { 1 Y2 };\n", \ "(0.5, 0.86) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Right 60 (1)\"\n", \ "\n", \ "quality 10\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0) { 0.5, 0, 1.0 };\n", \ "(0.5, 0.866) { 0.6, 0, 0.8 };\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(-0.125, 0.6495) { -0.5 X2, 0.75 X3 } { -0.5 Y2, 0.75 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(0.25, 0.433) { 0.5 X3 } { 0.5 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Vis A Vis (2)\"\n", \ "\n", \ "quality 2\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1);\n", \ "(0, 1);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 4) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.5) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y3, -0.25 Y4 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "(-0.5, 0.5) { -0.25 X2, -0.25 X3, 0.75 X4 } { -0.25 Y3, 0.75 Y4 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n", \ "(1, 1) { 1 X3 } { 1 Y3 };\n", \ "(0, 1) { 1 X4 } { 1 Y4 };\n", \ "(0, 0.5) { 0.5 X4 } { 0.5 Y4 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "orientations\n", \ "(1, 3, 4);\n", \ "(2, 3, 4);\n", \ "(1, 2, 3);\n", \ "(1, 2, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Triangle Vis A Vis (200)\"\n", \ "\n", \ "quality 200\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.2, 0.693) { 0.8 X2, 0.8 X3 } { 0.8 Y2, 0.8 Y3 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(-0.2, 0.693) { -0.6 X2, 0.8 X3 } { -0.6 Y2, 0.8 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.75, 0.433) { 0.5 X2, 0.5 X3 } { 0.5 Y2, 0.5 Y3 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(0.25, 0.433) { 0.5 X3 } { 0.5 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"2 h Vis A Vis (1)\"\n", \ "\n", \ "quality 3000\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1.732);\n", \ "(0, 1.732);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 0.25 X3, 0.25 X4 } { 0.25 Y2, 0.25 Y3, 0.25 Y4 };\n", \ "\n", \ "newlines\n", \ "(1, 5);\n", \ "(5, 4);\n", \ "(3, 5);\n", \ "(5, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { 1 Y2 };\n", \ "(1.5, 0.866) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y2, 0.75 Y3, -0.25 Y4 };\n", \ "(1, 1.732) { 1 X3 } { 1 Y3 };\n", \ "(0, 1.732) { 1 X4 } { 1 Y4 };\n", \ "(-0.5, 0.866) { 0.75 X4, -0.25 X2, -0.25 X3 } { 0.75 Y4, -0.25 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 5);\n", \ "(3, 4, 5);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ 0};


	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {1, 2, 3},
		  new int[] {0, 2, 4},
		  new int[] {0, 1, 5}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {1, 2, 5},
		  new int[] {1, 3, 6},
		  new int[] {1, 4, 7},
		  new int[] {2, 3, 8},
		  new int[] {2, 4, 9},
		  new int[] {3, 4, 10}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {2, 3, 4},
		  new int[] {1, 3, 5},
		  new int[] {1, 2, 6}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_reftab =
	  {
		  new int[] {1, 6, 5},
		  new int[] {2, 4, 6},
		  new int[] {3, 5, 4},
		  new int[] {6, 4, 5}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {1, 2, 5},
		  new int[] {2, 3, 6},
		  new int[] {3, 4, 7},
		  new int[] {1, 4, 8},
		  new int[] {5, 7, 9}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_reftab =
	  {
		  new int[] {1, 5, 9, 8},
		  new int[] {5, 2, 6, 9},
		  new int[] {8, 9, 7, 4},
		  new int[] {9, 6, 3, 7}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {1, 2, 5},
		  new int[] {1, 3, 6},
		  new int[] {1, 4, 7},
		  new int[] {2, 3, 8},
		  new int[] {2, 4, 9},
		  new int[] {3, 4, 10}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_reftab =
	  {
		  new int[] {1, 5, 6, 7},
		  new int[] {5, 2, 8, 9},
		  new int[] {6, 8, 3, 10},
		  new int[] {7, 9, 10, 4},
		  new int[] {5, 6, 7, 9},
		  new int[] {5, 6, 9, 8},
		  new int[] {6, 7, 9, 10},
		  new int[] {6, 8, 10, 9}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static bool[] Refine_reverse = {false, false, false, false, false, true, false, true};
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {1, 2, 9},
		  new int[] {3, 4, 10},
		  new int[] {4, 1, 11},
		  new int[] {2, 3, 12},
		  new int[] {5, 6, 13},
		  new int[] {7, 8, 14},
		  new int[] {8, 5, 15},
		  new int[] {6, 7, 16},
		  new int[] {1, 5, 17},
		  new int[] {2, 6, 18},
		  new int[] {3, 7, 19},
		  new int[] {4, 8, 20},
		  new int[] {2, 8, 21}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_fbetw =
	  {
		  new int[] {11, 12, 22},
		  new int[] {9, 10, 22},
		  new int[] {13, 14, 23},
		  new int[] {15, 16, 23},
		  new int[] {9, 13, 24},
		  new int[] {17, 18, 24},
		  new int[] {12, 16, 25},
		  new int[] {18, 19, 25},
		  new int[] {19, 20, 26},
		  new int[] {10, 14, 26},
		  new int[] {11, 15, 27},
		  new int[] {17, 20, 27}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_reftab =
	  {
		  new int[] {1, 9, 22, 11, 17, 24, 21, 27},
		  new int[] {9, 2, 12, 22, 24, 18, 25, 21},
		  new int[] {11, 22, 10, 4, 27, 21, 26, 20},
		  new int[] {22, 12, 3, 10, 21, 25, 19, 26},
		  new int[] {17, 24, 21, 27, 5, 13, 23, 15},
		  new int[] {24, 18, 25, 21, 13, 6, 16, 23},
		  new int[] {27, 21, 26, 20, 15, 23, 14, 8},
		  new int[] {21, 25, 19, 26, 23, 16, 7, 14}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_betw =
	  {
		  new int[] {3, 1, 7},
		  new int[] {1, 2, 8},
		  new int[] {3, 2, 9},
		  new int[] {6, 4, 10},
		  new int[] {4, 5, 11},
		  new int[] {6, 5, 12},
		  new int[] {1, 4, 13},
		  new int[] {3, 6, 14},
		  new int[] {2, 5, 15}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_fbetw =
	  {
		  new int[] {7, 10, 16},
		  new int[] {14, 13, 16},
		  new int[] {11, 8, 17},
		  new int[] {13, 15, 17},
		  new int[] {12, 9, 18},
		  new int[] {14, 15, 18}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] Refine_reftab =
	  {
		  new int[] {1, 8, 7, 13, 17, 16},
		  new int[] {7, 8, 9, 16, 17, 18},
		  new int[] {7, 9, 3, 16, 18, 14},
		  new int[] {8, 2, 9, 17, 15, 18},
		  new int[] {13, 17, 16, 4, 11, 10},
		  new int[] {16, 17, 18, 10, 11, 12},
		  new int[] {16, 18, 14, 10, 12, 6},
		  new int[] {17, 15, 18, 11, 5, 12}
	  };


	  internal static double CalcElementBadness(Array<Point2d> points, Element2d elem)
	  {
		// badness = sqrt(3) /36 * circumference^2 / area - 1 +
		//           h / li + li / h - 2

		Vec2d v12 = new Vec2d();
		Vec2d v13 = new Vec2d();
		Vec2d v23 = new Vec2d();
		double l12;
		double l13;
		double l23;
		double cir;
		double area;
		double c = ngsimd.GlobalMembers.sqrt(3.0) / 36;

		v12 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
		v13 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
		v23 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(2));

		l12 = v12.Length();
		l13 = v13.Length();
		l23 = v23.Length();

		cir = l12 + l13 + l23;
		area = 0.5 * (v12.X() * v13.Y() - v12.Y() * v13.X());
		if (area < 1e-6)
		{
		return 1e8;
		}

		if (testmode)
		{
		(*testout) << "l = " << l12 << " + " << l13 << " + " << l23 << " = " << cir << ", area = " << area << "\n";
		(*testout) << "shapeerr = " << 10 * (c * cir * cir / area - 1) << "\n" << "sizeerr = " << 1 / l12 + l12 + 1 / l13 + l13 + 1 / l23 + l23 - 6 << "\n";
		}

		return 10 * (c * cir * cir / area - 1) + 1 / l12 + l12 + 1 / l13 + l13 + 1 / l23 + l23 - 6;
	  }
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int ApplyRules_timer = NgProfiler.CreateTimer("meshing2::ApplyRules");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static bool ApplyRules_firsttime = true;
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int ApplyRules_timer1 = NgProfiler.CreateTimer("meshing2::ApplyRules 1");
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//extern double minother;
//C++ TO C# CONVERTER NOTE: 'extern' variable declarations are not required in C#:
	//extern double minwithoutother;


	  internal static double CalcElementBadness(Array<Point3d, PointIndex.BASE> points, Element elem)
	  {
	  double vol;
	  double l;
	  double l4;
	  double l5;
	  double l6;
	  if (elem.GetNP() != 4)
	  {
		  if (elem.GetNP() == 5)
		  {
		  double z = points[elem.PNum(5)].Z();
		  if (z > -1e-8)
		  {
			  return 1e8;
		  }
		  return (-1 / z) - z; //  - 2;
		  }
		  return 0;
	  }

	  Vec3d v1 = points[elem.PNum(2)] - points[elem.PNum(1)];
	  Vec3d v2 = points[elem.PNum(3)] - points[elem.PNum(1)];
	  Vec3d v3 = points[elem.PNum(4)] - points[elem.PNum(1)];

	  vol = - (Cross(new netgen.Vec3d(v1), new netgen.Vec3d(v2)) * v3);
	  l4 = Dist(points[elem.PNum(2)], points[elem.PNum(3)]);
	  l5 = Dist(points[elem.PNum(2)], points[elem.PNum(4)]);
	  l6 = Dist(points[elem.PNum(3)], points[elem.PNum(4)]);

	  l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

	  //  testout << "vol = " << vol << " l = " << l << endl;
	  if (vol < 1e-8)
	  {
		  return 1e10;
	  }
	  //  (*testout) << "l^3/vol = " << (l*l*l / vol) << endl;

	  double err = ngsimd.GlobalMembers.pow(l * l * l / vol, 1.0 / 3.0) / 12;
	  return err;
	  }
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	public static int ApplyRules_cnt = 0;
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	Array<int> ApplyRules_lpi(4);


	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_trig =
	  {
		  new int[] {1, 2, 3},
		  new int[] {0, 2, 4},
		  new int[] {0, 1, 5}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_quad6 =
	  {
		  new int[] {0, 1, 4},
		  new int[] {3, 2, 5}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_quad8 =
	  {
		  new int[] {0, 1, 4},
		  new int[] {3, 2, 5},
		  new int[] {0, 3, 6},
		  new int[] {1, 2, 7}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_trig =
	  {
		  new int[] {1, 2, 3},
		  new int[] {0, 2, 4},
		  new int[] {0, 1, 5}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_quad6 =
	  {
		  new int[] {0, 1, 4},
		  new int[] {3, 2, 5}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_quad8 =
	  {
		  new int[] {0, 1, 4},
		  new int[] {3, 2, 5},
		  new int[] {0, 3, 6},
		  new int[] {1, 2, 7}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_tet =
	  {
		  new int[] {0, 1, 4},
		  new int[] {0, 2, 5},
		  new int[] {0, 3, 6},
		  new int[] {1, 2, 7},
		  new int[] {1, 3, 8},
		  new int[] {2, 3, 9}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_prism =
	  {
		  new int[] {0, 2, 6},
		  new int[] {0, 1, 7},
		  new int[] {1, 2, 8},
		  new int[] {3, 5, 9},
		  new int[] {3, 4, 10},
		  new int[] {4, 5, 11}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_prism15 =
	  {
		  new int[] {0, 1, 6},
		  new int[] {0, 2, 7},
		  new int[] {1, 2, 8},
		  new int[] {0, 3, 9},
		  new int[] {1, 4, 10},
		  new int[] {2, 5, 11},
		  new int[] {3, 4, 12},
		  new int[] {3, 5, 13},
		  new int[] {4, 5, 14}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_pyramid =
	  {
		  new int[] {0, 1, 5},
		  new int[] {3, 2, 6},
		  new int[] {3, 0, 7},
		  new int[] {1, 2, 8},
		  new int[] {0, 4, 9},
		  new int[] {1, 4, 10},
		  new int[] {2, 4, 11},
		  new int[] {3, 4, 12}
	  };
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] MakeSecondOrder_betw_hex =
	  {
		  new int[] {0, 1, 8},
		  new int[] {2, 3, 9},
		  new int[] {3, 0, 10},
		  new int[] {1, 2, 11},
		  new int[] {4, 5, 12},
		  new int[] {6, 7, 13},
		  new int[] {7, 4, 14},
		  new int[] {5, 6, 15},
		  new int[] {0, 4, 16},
		  new int[] {1, 5, 17},
		  new int[] {2, 6, 18},
		  new int[] {3, 7, 19}
	  };


	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int[][] ValidateSecondOrder_betweentab =
	  {
		  new int[] {1, 2, 5},
		  new int[] {1, 3, 6},
		  new int[] {1, 4, 7},
		  new int[] {2, 3, 8},
		  new int[] {2, 4, 9},
		  new int[] {3, 4, 10}
	  };

	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int ImproveVolumeMesh_timer = NgProfiler.CreateTimer("MeshSmoothing 2D");

	  internal const double c_trig = 0.14433756; // sqrt(3.0) / 12
	  internal const double c_trig4 = 0.57735026; // sqrt(3.0) / 3


	  public static double CalcTriangleBadness(double x2, double x3, double y3, double metricweight, double h)
	  {
		// badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 
		// p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);

		double cir_2 = (x2 * x2 + x3 * x3 + y3 * y3 - x2 * x3);
		double area = x2 * y3;

		if (area <= 1e-24 * cir_2)
		{
		  return 1e10;
		}

		double badness = c_trig4 * cir_2 / area - 1;

		if (metricweight > 0)
		{
		// add:  metricweight * (area / h^2 + h^2 / area - 2)

		double areahh = area / (h * h);
		badness += metricweight * (areahh + 1 / areahh - 2);
		}
		return badness;
	  }

	  public static void CalcTriangleBadness(double x2, double x3, double y3, double metricweight, double h, ref double badness, ref double g1x, ref double g1y)
	  {
		// old: badness = sqrt(3.0) /36 * circumference^2 / area - 1 
		// badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 
		// p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);


		double cir_2 = 2 * (x2 * x2 + x3 * x3 + y3 * y3 - x2 * x3);
		double area = 0.5 * x2 * y3;

		if (area <= 1e-24 * cir_2)
		{
		g1x = 0;
		g1y = 0;
		badness = 1e10;
		return;
		}

		badness = c_trig * cir_2 / area - 1;

		double c1 = -2 * c_trig / area;
		double c2 = 0.5 * c_trig * cir_2 / (area * area);
		g1x = c1 * (x2 + x3) + c2 * y3;
		g1y = c1 * (y3) + c2 * (x2 - x3);

		if (metricweight > 0)
		{
		// area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
		// add:  metricweight * (area / h^2 + h^2 / area - 2)

		area = x2 * y3;
		double dareax1 = -y3;
		double dareay1 = x3 - x2;

		double areahh = area / (h * h);
		double fac = metricweight * (areahh - 1 / areahh) / area;

		badness += metricweight * (areahh + 1 / areahh - 2);
		g1x += fac * dareax1;
		g1y += fac * dareay1;
		}
	  }



	  public static double CalcTriangleBadness(Point < 3> p1, Point < 3> p2, Point < 3> p3, double metricweight, double h)
	  {
		// badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 

		Vec < 3> e12 = p2.functorMethod - p1.functorMethod;
		Vec < 3> e13 = p3.functorMethod - p1.functorMethod;
		Vec < 3> e23 = p3.functorMethod - p2.functorMethod;

		double cir_2 = e12.Length2() + e13.Length2() + e23.Length2();
		double area = 0.5 * Cross(e12, e13).Length();

		if (area <= 1e-24 * cir_2)
		{
		  return 1e10;
		}

		double badness = c_trig * cir_2 / area - 1;

		if (metricweight > 0)
		{
		// add:  metricweight * (area / h^2 + h^2 / area - 2)
			area *= 2; // optimum for (2 area) is h^2
			double areahh = area / (h * h);
		badness += metricweight * (areahh + 1 / areahh - 2);
		}

		return badness;
	  }

	  public static double CalcTriangleBadnessGrad(Point < 3> p1, Point < 3> p2, Point < 3> p3, ref Vec < 3> gradp1, double metricweight, double h)
	  {
		// badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1 

		Vec < 3> e12 = p2.functorMethod - p1.functorMethod;
		Vec < 3> e13 = p3.functorMethod - p1.functorMethod;
		Vec < 3> e23 = p3.functorMethod - p2.functorMethod;

		double cir_2 = e12.Length2() + e13.Length2() + e23.Length2();
		Vec < 3> varea = Cross(e12, e13);
		double area = 0.5 * varea.Length();

		Vec < 3> dcir_2 = (-2) * (e12 + e13);
		Vec < 3> darea = (0.25 / area) * Cross(p2.functorMethod - p3.functorMethod, varea);

		if (area <= 1e-24 * cir_2)
		{
			gradp1 = 0;
			return 1e10;
		}

		double badness = c_trig * cir_2 / area - 1;
		gradp1 = c_trig * (1.0 / area * dcir_2 - cir_2 / (area * area) * darea);

		if (metricweight > 0)
		{
		// add:  metricweight * (area / h^2 + h^2 / area - 2)
			area *= 2; // optimum for (2 area) is h^2

			double areahh = area / (h * h);
		badness += metricweight * (areahh + 1 / areahh - 2);

			gradp1.functorMethod += (2 * metricweight * (1 / (h * h) - (h * h) / (area * area))) * darea;
		}

		return badness;
	  }




	  public static double CalcTriangleBadness(Point < 3> p1, Point < 3> p2, Point < 3> p3, Vec < 3> n, double metricweight, double h)
	  {
		Vec < 3> v1 = p2.functorMethod - p1.functorMethod;
		Vec < 3> v2 = p3.functorMethod - p1.functorMethod;

		Vec < 3> e1 = v1;
		Vec < 3> e2 = v2;

		e1 -= (e1 * n.functorMethod) * n.functorMethod;
		e1 /= (e1.Length() + 1e-24);
		e2 = Cross(n.functorMethod, e1);

		return CalcTriangleBadness((e1 * v1), (e1 * v2), (e2 * v2), metricweight, h);
	  }







	  public static MeshOptimize2d dummy = new MeshOptimize2d();

	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int ImproveMesh_timer = NgProfiler.CreateTimer("MeshSmoothing 2D");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int ImproveMesh_timer1 = NgProfiler.CreateTimer("MeshSmoothing 2D start");
	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	  public static int ImproveMesh_timer2 = NgProfiler.CreateTimer("MeshSmoothing 2D - BFGS");





	//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	static Timer CalcTotalBad_t("CalcTotalBad");

	public static double CalcTotalBad(Mesh.T_POINTS points, Array<Element, 0, uint> elements, MeshingParameters mp)
	{
	//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
	//  static Timer t("CalcTotalBad");
	  RegionTimer reg = new RegionTimer(CalcTotalBad_t);

	  double sum = 0;
	  double elbad;

	  tets_in_qualclass.SetSize(20);
	  tets_in_qualclass = 0;

	  double teterrpow = mp.opterrpow;

	  for (int i = 1; i <= elements.Size(); i++)
	  {
		  elbad = ngsimd.GlobalMembers.pow(max2(CalcBad(points, elements.Get(i), 0, mp), 1e-10), 1 / teterrpow);

		  int qualclass = (int)20 / elbad + 1;
		  if (qualclass < 1)
		  {
			  qualclass = 1;
		  }
		  if (qualclass > 20)
		  {
			  qualclass = 20;
		  }
		  tets_in_qualclass.Elem(qualclass)++;

		  sum += elbad;
	  }
	  return sum;
	}

	public static int WrongOrientation(Mesh.T_POINTS points, Element el)
	{
	  Point3d p1 = points[el.PNum(1)];
	  Point3d p2 = points[el.PNum(2)];
	  Point3d p3 = points[el.PNum(3)];
	  Point3d p4 = points[el.PNum(4)];

	  Vec3d v1 = new Vec3d(p1, p2);
	  Vec3d v2 = new Vec3d(p1, p3);
	  Vec3d v3 = new Vec3d(p1, p4);
	  Vec3d n = new Vec3d();

	  Cross(v1, v2, n);
	  double vol = n * v3;

	  return (vol > 0);
	}
	public static string[] tetrules = {"tolfak 0.5\n", \ "\n", \ "rule \"Free Tetrahedron\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0);\n", \ "(0.5, 0.866, 0);\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.288, -0.816)\n", \ "	{ 0.333 X1, 0.333 X2, 0.333 X3 }\n", \ "	{ 0.333 Y1, 0.333 Y2, 0.333 Y3 } { };\n", \ "\n", \ "newfaces\n", \ "(4, 1, 2);\n", \ "(4, 2, 3);\n", \ "(4, 3, 1);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1.6 P4, -0.2 P1, -0.2 P2, -0.2 P3 };\n", \ "{ -0.5 P1, 0.5 P2, 0.5 P3, 0.5 P4 };\n", \ "{ 0.5 P1, -0.5 P2, 0.5 P3, 0.5 P4 };\n", \ "{ 0.5 P1, 0.5 P2, -0.5 P3, 0.5 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 60\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "flags c;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 0.866, 0) { 0.5 };\n", \ "(0.5, 0.288, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 4, 3);\n", \ "(4, 2, 3);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ -0.35 P1, 0.45 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, -0.35 P2, 0.45 P3, 0.45 P4 };\n", \ "{ -0.05 P1, -0.05 P2, 0.7 P3, 0.4 P4 };\n", \ "\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P2, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.65 P3, 0.35 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 60 with edge(1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "flags c;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.8 };\n", \ "(0.5, 0.866, 0) { 0.8 };\n", \ "(0.5, 0.288, -0.816) { 0.8 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "mapedges\n", \ "(3, 4);\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 4, 3);\n", \ "(4, 2, 3);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.4 P1, 0.4 P4, 0.4 P3, -0.2 P2 };\n", \ "{ 0.4 P2, 0.4 P4, 0.4 P3, -0.2 P1 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P1, 0.3333 P4, 0.3334 P3 };\n", \ "{ 0.3333 P2, 0.3333 P4, 0.3334 P3 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron Vis a Vis Point (1)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 0.866, 0) { 0.5 };\n", \ "(0.5, 0.288, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(4, 3, 1);\n", \ "(4, 2, 3);\n", \ "(4, 1, 2);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ -0.5 P1, 0.5 P2, 0.5 P3, 0.5 P4 };\n", \ "{ 0.5 P1, -0.5 P2, 0.5 P3, 0.5 P4 };\n", \ "{ 0.5 P1, 0.5 P2, -0.5 P3, 0.5 P4 };\n", \ "{ 0.8 P1, -0.1 P2, -0.1 P3, 0.4 P4 };\n", \ "{ -0.1 P1, 0.8 P2, -0.1 P3, 0.4 P4 };\n", \ "{ -0.1 P1, -0.1 P2, 0.8 P3, 0.4 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P2, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P2, 0.3334 P4 };\n", \ "{ 0.7 P1, 0.3 P4 };\n", \ "{ 0.7 P2, 0.3 P4 };\n", \ "{ 0.7 P3, 0.3 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron Vis a Vis Point with edge(1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 0.866, 0) { 0.5 };\n", \ "(0.5, 0.288, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "mapedges\n", \ "(1, 4);\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(4, 3, 1);\n", \ "(4, 2, 3);\n", \ "(4, 1, 2);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ -0.35 P1, 0.45 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, -0.35 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, 0.45 P2, -0.35 P3, 0.45 P4 };\n", \ "{ -0.05 P1, 0.7 P2, -0.05 P3, 0.4 P4 };\n", \ "{ -0.05 P1, -0.05 P2, 0.7 P3, 0.4 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P2, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P2, 0.3334 P4 };\n", \ "{ 0.65 P2, 0.35 P4 };\n", \ "{ 0.65 P3, 0.35 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron Vis a Vis Point with 2 edges (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 0.866, 0) { 0.5 };\n", \ "(0.5, 0.288, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "mapedges\n", \ "(1, 4);\n", \ "(2, 4);\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(4, 3, 1);\n", \ "(4, 2, 3);\n", \ "(4, 1, 2);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ -0.35 P1, 0.45 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, -0.35 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, 0.45 P2, -0.35 P3, 0.45 P4 };\n", \ "{ -0.05 P1, -0.05 P2, 0.7 P3, 0.4 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P2, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P2, 0.3334 P4 };\n", \ "{ 0.65 P3, 0.35 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron Vis a Vis Point with 3 edges (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 0.866, 0) { 0.5 };\n", \ "(0.5, 0.288, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "mapedges\n", \ "(1, 4);\n", \ "(2, 4);\n", \ "(3, 4);\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(4, 3, 1);\n", \ "(4, 2, 3);\n", \ "(4, 1, 2);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ -0.35 P1, 0.45 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, -0.35 P2, 0.45 P3, 0.45 P4 };\n", \ "{ 0.45 P1, 0.45 P2, -0.35 P3, 0.45 P4 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P2, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P3, 0.3334 P4 };\n", \ "{ 0.3333 P1, 0.3333 P2, 0.3334 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron Vis a Vis Triangle (1)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 0.866, 0) { 0.5 };\n", \ "(0, 0, -0.816) { 0.5 };\n", \ "(1, 0, -0.816) { 0.5 };\n", \ "(0.5, 0.866, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(4, 6, 5) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 2, 4);\n", \ "(2, 5, 4);\n", \ "(2, 3, 6);\n", \ "(2, 6, 5);\n", \ "(3, 1, 4);\n", \ "(3, 4, 6);\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "(4, 2, 3, 6);\n", \ "(4, 2, 6, 5);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ -0.2 P1, 0.35 P2, 0.35 P3, -0.2 P4, 0.35 P5, 0.35 P6 };\n", \ "{ 0.35 P1, -0.2 P2, 0.35 P3, 0.35 P4, -0.2 P5, 0.35 P6 };\n", \ "{ 0.35 P1, 0.35 P2, -0.2 P3, 0.35 P4, 0.35 P5, -0.2 P6 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Octaeder 1\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.95 };\n", \ "(0.5, 0.866, 0) { 0.95 };\n", \ "(0.5, -0.288, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "newpoints\n", \ "(1, 0.578, -0.816) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "(0, 0.578, -0.816) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "\n", \ "newfaces\n", \ "(2, 3, 5);\n", \ "(3, 1, 6);\n", \ "(3, 6, 5);\n", \ "(2, 5, 4);\n", \ "(1, 4, 6);\n", \ "(4, 5, 6);\n", \ "\n", \ "elements\n", \ "(3, 4, 1, 2);\n", \ "(3, 4, 2, 5);\n", \ "(3, 4, 5, 6);\n", \ "(3, 4, 6, 1);\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, -0.288, -0.816) { 1 X4 } { 1 Y4 } { 1 Z4 };\n", \ "(-0.5, 1, -1.5) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 } { 1 Z4 };\n", \ "( 1.5, 1, -1.5) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 } { 1 Z4 };\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Octaeder 2\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.95 };\n", \ "(0.5, 0.866, 0) { 0.95 };\n", \ "(0.5, -0.288, -0.816) { 0.5 };\n", \ "(1, 0.578, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0, 0.578, -0.816) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "\n", \ "newfaces\n", \ "(2, 3, 5);\n", \ "(3, 1, 6);\n", \ "(3, 6, 5);\n", \ "(2, 5, 4);\n", \ "(1, 4, 6);\n", \ "(4, 5, 6);\n", \ "\n", \ "elements\n", \ "(3, 4, 1, 2);\n", \ "(3, 4, 2, 5);\n", \ "(3, 4, 5, 6);\n", \ "(3, 4, 6, 1);\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, -0.288, -0.816) { 1 X4 } { 1 Y4 } { 1 Z4 };\n", \ "(1, 0.578, -0.816) { 1 X5 } { 1 Y5 } { 1 Z5 };\n", \ "\n", \ "(0.9, 0.097, -0.544) { 0.333 X2, 0.333 X4, 0.333 X5 }\n", \ "                     { 0.333 Y2, 0.333 Y4, 0.333 Y5 }\n", \ "                     { 0.333 Z2, 0.333 Z4, 0.333 Z5 };\n", \ "(0.9, 0.481, -0.272) { 0.333 X2, 0.333 X3, 0.333 X5 }\n", \ "                     { 0.333 Y2, 0.333 Y3, 0.333 Y5 }\n", \ "                     { 0.333 Z2, 0.333 Z3, 0.333 Z5 };\n", \ "\n", \ "(-0.5, 1, -1.5) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 } { 0.5 Z4, 0.5 Z5 };\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"Octaeder 2a\"\n", \ "\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.95 };\n", \ "(0.5, 0.866, 0) { 0.95 };\n", \ "(0.5, -0.288, -0.816) { 0.5 };\n", \ "(1, 0.578, -0.816) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(3, 2, 5) del;\n", \ "\n", \ "newpoints\n", \ "(0, 0.578, -0.816)\n", \ "	{ -1 X2, 1 X3, 1 X4 }\n", \ "	{ -1 Y2, 1 Y3, 1 Y4 }\n", \ "	{ -1 Z2, 1 Z3, 1 Z4 };\n", \ "\n", \ "newfaces\n", \ "(1, 2, 4);\n", \ "(3, 1, 6);\n", \ "(3, 6, 5);\n", \ "(2, 5, 4);\n", \ "(1, 4, 6);\n", \ "(4, 5, 6);\n", \ "\n", \ "elements\n", \ "(3, 4, 1, 2);\n", \ "(3, 4, 2, 5);\n", \ "(3, 4, 5, 6);\n", \ "(3, 4, 6, 1);\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, -0.288, -0.816) { 1 X4 } { 1 Y4 } { 1 Z4 };\n", \ "(1, 0.578, -0.816) { 1 X5 } { 1 Y5 } { 1 Z5 };\n", \ "\n", \ "(0.9, 0.097, -0.544) { 0.333 X2, 0.333 X4, 0.333 X5 }\n", \ "                     { 0.333 Y2, 0.333 Y4, 0.333 Y5 }\n", \ "                     { 0.333 Z2, 0.333 Z4, 0.333 Z5 };\n", \ "\n", \ "(0.5, -0.097, -0.272) { 0.333 X2, 0.333 X4, 0.333 X1 }\n", \ "                     { 0.333 Y2, 0.333 Y4, 0.333 Y1 }\n", \ "                     { 0.333 Z2, 0.333 Z4, 0.333 Z1 };\n", \ "\n", \ "(-0.5, 1, -1.5) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 } { 0.5 Z4, 0.5 Z5 };\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Pyramid 1\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(0.5, 0.866, 0) { 1 };\n", \ "(0.5, -0.288, -0.816) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "newpoints\n", \ "(1, 0.578, -0.816) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "\n", \ "newfaces\n", \ "(1, 4, 3);\n", \ "(2, 3, 5);\n", \ "(2, 5, 4);\n", \ "(4, 5, 3);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "(4, 2, 3, 5);\n", \ "\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, -0.288, -0.816) { 1 X4 } { 1 Y4 } { 1 Z4 };\n", \ "(0, 1, -1) { 0.5 X3, 0.5 X4 } { 1 Y3 } { 1 Z4 };\n", \ "(1.5, 1, -1.5) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 2 times 60\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.3 };\n", \ "(0.5, 0.866, 0) { 0.3 };\n", \ "(0.5, 0.288, -0.816) { 0.3 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(2, 4, 3) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "(1, 4, 3);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.4 P1, 0.4 P4, 0.4 P3, -0.2 P2 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 0.3333 P1, 0.3333 P3, 0.3334 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Fill Tetrahedron (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.2 };\n", \ "(0.5, 0.866, 0) { 0.2 };\n", \ "(0.5, 0.288, -0.816) { 0.2 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(2, 4, 3) del;\n", \ "(3, 4, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newfaces\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 120 (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(0.5, 0.866, 0) { 1 };\n", \ "(0.5, -0.674, -0.544) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.288, -0.816)\n", \ "	{ -0.5 X1, -0.5 X2, 1 X3, 1 X4 }\n", \ "	{ -0.5 Y1, -0.5 Y2, 1 Y3, 1 Y4}\n", \ "	{ -0.5 Z1, -0.5 Z2, 1 Z3, 1 Z4};\n", \ "\n", \ "newfaces\n", \ "(1, 5, 3);\n", \ "(3, 5, 2);\n", \ "(1, 4, 5);\n", \ "(2, 5, 4);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 5);\n", \ "(1, 4, 2, 5);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.3 P5, -0.3 P1 };\n", \ "{ 1.3 P5, -0.3 P2 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P5 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 2 times 120 (1)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(0.5, 0.866, 0) { 1 };\n", \ "(0.5, -0.674, -0.544) { 0.8 };\n", \ "(1.334, 0.77, -0.544) { 0.8 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(3, 2, 5) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.288, -0.816) { 0.25 X1, -0.5 X2, 0.25 X3, 0.5 X4, 0.5 X5 }\n", \ "		 { 0.25 Y1, -0.5 Y2, 0.25 Y3, 0.5 Y4, 0.5 Y5 }\n", \ "		 { 0.25 Z1, -0.5 Z2, 0.25 Z3, 0.5 Z4, 0.5 Z5 };\n", \ "\n", \ "newfaces\n", \ "(6, 3, 1);\n", \ "(6, 1, 4);\n", \ "(6, 4, 2);\n", \ "(6, 2, 5);\n", \ "(6, 5, 3);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 6);\n", \ "(1, 4, 2, 6);\n", \ "(2, 5, 3, 6);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1.4 P6, -0.4 P2 };\n", \ "{ 1.4 P6, -0.4 P1 };\n", \ "{ 1.4 P6, -0.4 P3 };\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"four Tetrahedron non convex (4)\"\n", \ "\n", \ "quality 4\n", \ "flags l;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.1 };\n", \ "(0.5, 1, 0) { 0.1 };\n", \ "(0.5, 0, -1) { 0.1 };\n", \ "(0.5, 0.3, -0.3) { 0.1 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(1, 5, 4) del;\n", \ "(1, 3, 5) del;\n", \ "\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.1, -0.1)\n", \ "	 { 0.333 X1, 0.333 X2, 0.334 X5 }\n", \ "	 { 0.333 Y1, 0.333 Y2, 0.334 Y5 }\n", \ "	 { 0.333 Z1, 0.333 Z2, 0.334 Z5 };\n", \ "\n", \ "newfaces\n", \ "(6, 2, 3) del;\n", \ "(6, 4, 2) del;\n", \ "(6, 5, 4) del;\n", \ "(6, 3, 5) del;\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 6);\n", \ "(1, 4, 2, 6);\n", \ "(1, 5, 4, 6);\n", \ "(1, 3, 5, 6);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1.5 P6, -0.5 P1 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "\n", \ "\n", \ "freeset\n", \ "1 6 2 3;\n", \ "\n", \ "freeset\n", \ "1 6 3 5;\n", \ "\n", \ "freeset\n", \ "1 6 5 4;\n", \ "\n", \ "freeset\n", \ "1 6 4 2;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"five Tetrahedron non convex (4)\"\n", \ "\n", \ "quality 4\n", \ "flags l;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0, 0.8, -0.2) { 0.5 };\n", \ "(0, 0.2, -0.8) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 3, 4) del;\n", \ "(1, 4, 5) del;\n", \ "(1, 5, 6) del;\n", \ "(1, 6, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.1, 0.1, -0.1)\n", \ "	 { 0.75 X1, 0.05 X2, 0.05 X3, 0.05 X4, 0.05 X5, 0.05 X6 }\n", \ "	 { 0.75 Y1, 0.05 Y2, 0.05 Y3, 0.05 Y4, 0.05 Y5, 0.05 Y6 }\n", \ "	 { 0.75 Z1, 0.05 Z2, 0.05 Z3, 0.05 Z4, 0.05 Z5, 0.05 Z6 };\n", \ "\n", \ "newfaces\n", \ "(7, 2, 3);\n", \ "(7, 3, 4);\n", \ "(7, 4, 5);\n", \ "(7, 5, 6);\n", \ "(7, 6, 2);\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 7);\n", \ "(1, 3, 4, 7);\n", \ "(1, 4, 5, 7);\n", \ "(1, 5, 6, 7);\n", \ "(1, 6, 2, 7);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1.5 P7, -0.5 P1 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "{ 1 P7 };\n", \ "\n", \ "\n", \ "\n", \ "freeset\n", \ "1 7 2 3;\n", \ "\n", \ "freeset\n", \ "1 7 3 4;\n", \ "\n", \ "freeset\n", \ "1 7 4 5;\n", \ "\n", \ "freeset\n", \ "1 7 5 6;\n", \ "\n", \ "freeset\n", \ "1 7 6 2;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"four Tetrahedron non convex (6)\"\n", \ "\n", \ "quality 6\n", \ "flags l;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "(0.5, 0.3, -0.3) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(1, 5, 4) del;\n", \ "(1, 3, 5) del;\n", \ "\n", \ "\n", \ "newpoints\n", \ "(0.095, 0.003, -0.003)\n", \ "	 { 0.9 X1, 0.09 X2, 0.01 X5 }\n", \ "	 { 0.9 Y1, 0.09 Y2, 0.01 Y5 }\n", \ "	 { 0.9 Z1, 0.09 Z2, 0.01 Z5 };\n", \ "\n", \ "newfaces\n", \ "(6, 2, 3) del;\n", \ "(6, 4, 2) del;\n", \ "(6, 5, 4) del;\n", \ "(6, 3, 5) del;\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 6);\n", \ "(1, 4, 2, 6);\n", \ "(1, 5, 4, 6);\n", \ "(1, 3, 5, 6);\n", \ "\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1.499 P6, -0.5 P1, 0.001 P2 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "\n", \ "\n", \ "freeset\n", \ "1 6 2 3;\n", \ "\n", \ "freeset\n", \ "1 6 3 5;\n", \ "\n", \ "freeset\n", \ "1 6 5 4;\n", \ "\n", \ "freeset\n", \ "1 6 4 2;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"four Tetrahedron non convex (6)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "(0.5, 0.4, -0.4) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(4, 5, 2) del;\n", \ "(5, 3, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.925, 0.02, -0.02)\n", \ "	 { 0.05 X1, 0.9 X2, 0.05 X5 }\n", \ "	 { 0.05 Y1, 0.9 Y2, 0.05 Y5 }\n", \ "	 { 0.05 Z1, 0.9 Z2, 0.05 Z5 };\n", \ "\n", \ "newfaces\n", \ "(3, 1, 6);\n", \ "(1, 4, 6);\n", \ "(4, 5, 6);\n", \ "(5, 3, 6);\n", \ "\n", \ "elements\n", \ "(3, 1, 2, 6);\n", \ "(1, 4, 2, 6);\n", \ "(4, 5, 2, 6);\n", \ "(5, 3, 2, 6);\n", \ "\n", \ "orientations\n", \ "(3, 1, 2, 5);\n", \ "(1, 4, 2, 5);\n", \ "(2, 4, 5, 1);\n", \ "(3, 2, 5, 1);\n", \ "(5, 4, 2, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1.5 P6, -0.5 P2 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "{ 1 P6 };\n", \ "\n", \ "freeset\n", \ "3 1 2 6;\n", \ "\n", \ "freeset\n", \ "1 4 2 6;\n", \ "\n", \ "freeset\n", \ "4 5 2 6;\n", \ "\n", \ "freeset\n", \ "5 3 2 6;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"three Tetrahedron non convex (4)\"\n", \ "\n", \ "quality 4\n", \ "flags l;\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(1, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.25, -0.25)\n", \ "	 { 0.25 X1, 0.25 X2, 0.25 X3, 0.25 X4 }\n", \ "	 { 0.25 Y1, 0.25 Y2, 0.25 Y3, 0.25 Y4 }\n", \ "	 { 0.25 Z1, 0.25 Z2, 0.25 Z3, 0.25 Z4 };\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3);\n", \ "(5, 4, 2);\n", \ "(5, 3, 4);\n", \ "\n", \ "elements\n", \ "(2, 3, 1, 5);\n", \ "(3, 4, 1, 5);\n", \ "(4, 2, 1, 5;\n", \ "\n", \ "orientations\n", \ "(1, 2, 4, 3);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.5 P5, -0.5 P1 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "freeset\n", \ "1 4 2 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"three Tetrahedron non convex (6)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(1, 3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.2, 0.1, -0.1)\n", \ "	 { 0.7 X1, 0.1 X2, 0.1 X3, 0.1 X4 }\n", \ "	 { 0.7 Y1, 0.1 Y2, 0.1 Y3, 0.1 Y4 }\n", \ "	 { 0.7 Z1, 0.1 Z2, 0.1 Z3, 0.1 Z4 };\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3);\n", \ "(5, 4, 2);\n", \ "(5, 3, 4);\n", \ "\n", \ "elements\n", \ "(2, 3, 1, 5);\n", \ "(3, 4, 1, 5);\n", \ "(4, 2, 1, 5;\n", \ "\n", \ "orientations\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.5 P5, -0.5 P1 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 3 4 5;\n", \ "\n", \ "freeset\n", \ "1 4 2 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"four Tetrahedron non convex (6)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "(0.5, 0.4, -0.4) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "(4, 5, 2) del;\n", \ "(5, 3, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.7, 0.08, -0.08) { 0.6 X2, 0.2 X5 } { 0.2 Y5 } { 0.2 Z5 };\n", \ "\n", \ "newfaces\n", \ "(3, 1, 6);\n", \ "(1, 4, 6);\n", \ "(4, 5, 6);\n", \ "(5, 3, 6);\n", \ "\n", \ "elements\n", \ "(3, 1, 2, 6);\n", \ "(1, 4, 2, 6);\n", \ "(4, 5, 2, 6);\n", \ "(5, 3, 2, 6);\n", \ "\n", \ "\n", \ "orientations\n", \ "(3, 1, 2, 5);\n", \ "(5, 1, 2, 4);\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 1, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, 0, -1) { 1 X4 } { 1 Y4 } { 1 Z4 };\n", \ "(0.5, 0.4, -0.4) { 1 X5 } { 1 Y5 } { 1 Z5 };\n", \ "(0.55, 0.12, -0.12) { 0.4 X2, 0.3 X5 } { 0.3 Y5 } { 0.3 Z5 };\n", \ "\n", \ "freeset\n", \ "3 1 2 6;\n", \ "\n", \ "freeset\n", \ "1 4 2 6;\n", \ "\n", \ "freeset\n", \ "4 5 2 6;\n", \ "\n", \ "freeset\n", \ "5 3 2 6;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 2 in 60 (12)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 0.5 };\n", \ "(0.5, 1, 0) { 0.5 };\n", \ "(0.5, 0, -1) { 0.5 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.1, -0.1)\n", \ "	{ 0.4 X1, 0.4 X2, 0.1 X3, 0.1 X4 }\n", \ "	{ 0.4 Y1, 0.4 Y2, 0.1 Y3, 0.1 Y4 }\n", \ "	{ 0.4 Z1, 0.4 Z2, 0.1 Z3, 0.1 Z4 };\n", \ "\n", \ "newfaces\n", \ "(5, 2, 3);\n", \ "(5, 3, 1);\n", \ "(5, 4, 2);\n", \ "(5, 1, 4);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 5);\n", \ "(1, 2, 5, 4);\n", \ "\n", \ "freezone2\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1.5 P5, -0.25 P1, -0.25 P2 };\n", \ "\n", \ "freezonelimit\n", \ "{ 1 P1 };\n", \ "{ 1 P2 };\n", \ "{ 1 P3 };\n", \ "{ 1 P4 };\n", \ "{ 1 P5 };\n", \ "\n", \ "freeset\n", \ "1 2 3 5;\n", \ "\n", \ "freeset\n", \ "1 2 4 5;\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Tetrahedron 120, but more than 180 (13)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 };\n", \ "(0.5, 0.866, 0) { 1 };\n", \ "(0.5, -0.866, 0) { 1 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "(1, 4, 2);\n", \ "\n", \ "newpoints\n", \ "(0.5, 0, -0.3) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "\n", \ "newfaces\n", \ "(1, 5, 3);\n", \ "(3, 5, 2);\n", \ "(2, 5, 1);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 5);\n", \ "\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, -0.1, -0.4)  { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4} { 0.5 Z3, 0.5 Z4 };\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Free Tetrahedron (14)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1.0 };\n", \ "(0.5, 0.866, 0) { 1.0 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.288, -0.2) { 0.333 X2, 0.333 X3 } { 0.333 Y3 } { };\n", \ "\n", \ "newfaces\n", \ "(4, 1, 2);\n", \ "(4, 2, 3);\n", \ "(4, 3, 1);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, 0.288, -0.25) { 0.333 X2, 0.333 X3 } { 0.333 Y3 } { };\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Free Tetrahedron (15)\"\n", \ "\n", \ "quality 100\n", \ "\n", \ "mappoints\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1.0 };\n", \ "(0.5, 0.866, 0) { 1.0 };\n", \ "\n", \ "mapfaces\n", \ "(1, 2, 3) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.288, -0.1) { 0.333 X2, 0.333 X3 } { 0.333 Y3 } { };\n", \ "\n", \ "newfaces\n", \ "(4, 1, 2);\n", \ "(4, 2, 3);\n", \ "(4, 3, 1);\n", \ "\n", \ "elements\n", \ "(1, 2, 3, 4);\n", \ "\n", \ "freezone\n", \ "(0, 0, 0);\n", \ "(1, 0, 0) { 1 X2 } { } { };\n", \ "(0.5, 0.866, 0) { 1 X3 } { 1 Y3 } { };\n", \ "(0.5, 0.288, -0.15) { 0.333 X2, 0.333 X3 } { 0.333 Y3 } { };\n", \ "\n", \ "endrule\n", null};
	  public static void QuickSortRec<T>(FlatArray<T> data, int left, int right)
	  {
		int i = left;
		int j = right;
		T midval = data[(left + right) / 2];

		do
		{
		while (data[i] < midval)
		{
			i++;
		}
		while (midval < data[j])
		{
			j--;
		}

		if (i <= j)
		{
			Swap(ref data[i], ref data[j]);
			i++;
			j--;
		}
		} while (i <= j);
		if (left < j)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: QuickSortRec(data, left, j);
			QuickSortRec(new FlatArray<T>(data), left, j);
		}
		if (i < right)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: QuickSortRec(data, i, right);
			QuickSortRec(new FlatArray<T>(data), i, right);
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class T>
	  public static void QuickSort<T>(FlatArray<T> data)
	  {
		if (data.Size() > 1)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: QuickSortRec(data, 0, data.Size()-1);
		  QuickSortRec(new FlatArray<T>(data), 0, (int)(data.Size() - 1));
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename FUNC>
	  public static void LoopOverEdges<FUNC>(Mesh mesh, MeshTopology top, PointIndex v, FUNC func)
	  {
		foreach (ElementIndex elnr in top.GetVertexElements(new netgen.PointIndex(v)))
		{
			Element el = mesh[elnr];

			int neledges = MeshTopology.GetNEdges(el.GetType());
			ELEMENT_EDGE[] eledges = MeshTopology.GetEdges0(el.GetType());

			for (int k = 0; k < neledges; k++)
			{
				INDEX_2 edge = new INDEX_2(el[eledges[k][0]], el[eledges[k][1]]);
				// edge.Sort();
				int edgedir = (edge.I1() > edge.I2());
				if (edgedir != 0)
				{
					swap(edge.I1(), edge.I2());
				}
				if (edge.I1() != v)
				{
					continue;
				}

				func(edge, elnr, k, 3, edgedir);
			}
		}

		foreach (SurfaceElementIndex elnr in top.GetVertexSurfaceElements(new netgen.PointIndex(v)))
		{
			Element2d el = mesh[elnr];

			int neledges = MeshTopology.GetNEdges(el.GetType());
			ELEMENT_EDGE[] eledges = MeshTopology.GetEdges0(el.GetType());

			for (int k = 0; k < neledges; k++)
			{
				INDEX_2 edge = new INDEX_2(el[eledges[k][0]], el[eledges[k][1]]);
				// edge.Sort();
				int edgedir = (edge.I1() > edge.I2());
				if (edgedir != 0)
				{
					swap(edge.I1(), edge.I2());
				}

				if (edge.I1() != v)
				{
					continue;
				}

				func(edge, elnr, k, 2, edgedir);
			}
		}

		foreach (SegmentIndex elnr in top.GetVertexSegments(new netgen.PointIndex(v)))
		{
			Segment el = mesh[elnr];
			INDEX_2 edge = new INDEX_2(el[0], el[1]);
			int edgedir = (edge.I1() > edge.I2());
			if (edgedir != 0)
			{
				swap(edge.I1(), edge.I2());
			}

			edge.Sort();
			if (edge.I1() != v)
			{
				continue;
			}

			func(edge, elnr, 0, 1, edgedir);
		}
	  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename FUNC>
	  public static void LoopOverFaces<FUNC>(Mesh mesh, MeshTopology top, PointIndex v, FUNC func)
	  {
		foreach (ElementIndex elnr in top.GetVertexElements(new netgen.PointIndex(v)))
		{
			Element el = mesh[elnr];

			int nelfaces = MeshTopology.GetNFaces(el.GetType());
			ELEMENT_FACE[] elfaces = MeshTopology.GetFaces0(el.GetType());

			for (int j = 0; j < nelfaces; j++)
			{
			  if (elfaces[j][3] < 0)

			  { // triangle
				  INDEX_4 face = new INDEX_4(el[elfaces[j][0]], el[elfaces[j][1]], el[elfaces[j][2]], 0);

				  int facedir = 0;
				  if (face.I1() > face.I2())
				  {
						swap(face.I1(), face.I2());
						facedir += 1;
				  }
				  if (face.I2() > face.I3())
				  {
						swap(face.I2(), face.I3());
						facedir += 2;
				  }
				  if (face.I1() > face.I2())
				  {
						swap(face.I1(), face.I2());
						facedir += 4;
				  }

				  if (face.I1() != v)
				  {
					  continue;
				  }

				  func(face, elnr, j, true, facedir);
			  }
			/*
			      if (pass == 1)
			        {
			          if (!vert2face.Used (face))
			            {
			              nfa++;
			              vert2face.Set (face, nfa);
			              INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
			              face2vert.Append (hface);
			            }
			        }
			      else
			        {
			          int facenum = vert2face.Get(face);
			          faces[elnr][j].fnr = facenum-1;
			          faces[elnr][j].forient = facedir;
			        }
			*/
			  else
			  {
				  // quad
				  // int facenum;
				  INDEX_4 face4 = new INDEX_4(el[elfaces[j][0]], el[elfaces[j][1]], el[elfaces[j][2]], el[elfaces[j][3]]);

				  int facedir = 0;
				  if (min2(face4.I1(), face4.I2()) > min2(face4.I4(), face4.I3()))
				  { // z - flip
					  facedir += 1;
					  swap(face4.I1(), face4.I4());
					  swap(face4.I2(), face4.I3());
				  }
				  if (min2(face4.I1(), face4.I4()) > min2(face4.I2(), face4.I3()))
				  { // x - flip
					  facedir += 2;
					  swap(face4.I1(), face4.I2());
					  swap(face4.I3(), face4.I4());
				  }
				  if (face4.I2() > face4.I4())
				  { // diagonal flip
					  facedir += 4;
					  swap(face4.I2(), face4.I4());
				  }

				  if (face4.I1() != v)
				  {
					  continue;
				  }

				  func(face4, elnr, j, true, facedir);
					/*
				  INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
	
	
				  if (vert2face.Used (face))
					{
					  facenum = vert2face.Get(face);
					}
				  else
					{
					  if (pass == 2) cout << "hier in pass 2" << endl;
					  nfa++;
					  vert2face.Set (face, nfa);
					  facenum = nfa;
					  
					  INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
					  face2vert.Append (hface);
					}
	
				  faces[elnr][j].fnr = facenum-1;
					          faces[elnr][j].forient = facedir;
				}
					*/
			  }
			}
		}


		foreach (SurfaceElementIndex elnr in top.GetVertexSurfaceElements(new netgen.PointIndex(v)))
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: const Element2d & el = mesh.SurfaceElement(elnr);
			Element2d el = mesh.SurfaceElement(new netgen.SurfaceElementIndex(elnr));

			ELEMENT_FACE[] elfaces = MeshTopology.GetFaces1(el.GetType());

			if (elfaces[0][3] == 0)

			{ // triangle

				// int facenum;
				int facedir;

				INDEX_4 face = new INDEX_4(el.PNum(elfaces[0][0]), el.PNum(elfaces[0][1]), el.PNum(elfaces[0][2]), 0);

				facedir = 0;
				if (face.I1() > face.I2())
				{
					swap(face.I1(), face.I2());
					facedir += 1;
				}
				if (face.I2() > face.I3())
				{
					swap(face.I2(), face.I3());
					facedir += 2;
				}
				if (face.I1() > face.I2())
				{
					swap(face.I1(), face.I2());
					facedir += 4;
				}

				if (face.I1() != v)
				{
					continue;
				}

				func(face, elnr, 0, false, facedir);
				/*
				  if (vert2face.Used (face))
				  facenum = vert2face.Get(face);
				  else
				  {
				  nfa++;
				  vert2face.Set (face, nfa);
				  facenum = nfa;
				  
				  INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
				  face2vert.Append (hface);
				  }
				    
				  surffaces[elnr].fnr = facenum-1;
				  surffaces[elnr].forient = facedir;
				*/
			}

			else

			{
				// quad
				// int facenum;
				int facedir;

				INDEX_4 face4 = new INDEX_4(el.PNum(elfaces[0][0]), el.PNum(elfaces[0][1]), el.PNum(elfaces[0][2]), el.PNum(elfaces[0][3]));

				facedir = 0;
				if (min2(face4.I1(), face4.I2()) > min2(face4.I4(), face4.I3()))
				{ // z - orientation
					facedir += 1;
					swap(face4.I1(), face4.I4());
					swap(face4.I2(), face4.I3());
				}
				if (min2(face4.I1(), face4.I4()) > min2(face4.I2(), face4.I3()))
				{ // x - orientation
					facedir += 2;
					swap(face4.I1(), face4.I2());
					swap(face4.I3(), face4.I4());
				}
				if (face4.I2() > face4.I4())
				{
					facedir += 4;
					swap(face4.I2(), face4.I4());
				}

				if (face4.I1() != v)
				{
					continue;
				}
				func(face4, elnr, 0, false, facedir);
				/*
				  INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
	
				  if (vert2face.Used (face))
				  facenum = vert2face.Get(face);
				  else
				  {
				  nfa++;
				  vert2face.Set (face, nfa);
				  facenum = nfa;
				        
				  INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
				  face2vert.Append (hface);
				  }
				    
				  surffaces[elnr].fnr = facenum-1;
				  surffaces[elnr].forient = facedir;
				  }
				*/
			}

		}
	  }
	public static string[] triarules = {"rule \"Free Triangle (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0) { 1.0, 0, 1.0 };\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 0.5 X2 } { };\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.7) { 0.5 X2 } { };\n", \ "(0.5, 1.5) { 0.5 X2 } { };\n", \ "(-0.5, 0.7) { 0.5 X2 } { };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.866) { 0.5 X2 } { };\n", \ "(0.5, 0.866) { 0.5 X2 } { };\n", \ "(0.5, 0.866) { 0.5 X2 } { };\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"Free Triangle (5)\"\n", \ "\n", \ "quality 5\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0) { 1.0, 0, 1.0 };\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.5) { 0.5 X2 } { };\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.7) { 1 X2 } { };\n", \ "(0, 0.7) { } { };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.5) { 0.5 X2 } { };\n", \ "(0.5, 0.5) { 0.5 X2 } { };\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Free Triangle (10)\"\n", \ "\n", \ "quality 10\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0) { 1.0, 0, 1.0 };\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.3) { 0.5 X2 } { };\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.5) { 1 X2 } { };\n", \ "(0, 0.5) { } { };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.3) { 0.5 X2 } { };\n", \ "(0.5, 0.3) { 0.5 X2 } { };\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Free Triangle (20)\"\n", \ "\n", \ "quality 20\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0) { 1.0, 0, 1.0 };\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.1) { 0.5 X2 } { };\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1, 0.2) { 1 X2 } { };\n", \ "(0, 0.2) { } { };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.1) { 0.5 X2 } { };\n", \ "(0.5, 0.1) { 0.5 X2 } { };\n", \ "\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Right 60 (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0) { 0.5, 0, 1.0 };\n", \ "(0.5, 0.866) { 0.6, 0, 0.8 };\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(-0.125, 0.6495) { -0.5 X2, 0.75 X3 } { -0.5 Y2, 0.75 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(0.25, 0.433) { 0.5 X3 } { 0.5 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left 60 (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.125, 0.6495) { 0.75 X2, 0.75 X3 } { 0.75 Y3 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.75, 0.433) { 0.5 X2, 0.5 X3 } { 0.5 Y2, 0.5 Y3 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Right 120 (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 1 X3, -1 X2 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(1, 4);\n", \ "(4, 3);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(1, 1.732) { -2 X2, 2 X3 } { 2 Y3 };\n", \ "(0, 1.732) { -3 X2, 2 X3 } { 2 Y3 };\n", \ "(-0.5, 0.866) { -2 X2, 1 X3 } {1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4);\n", \ "(2, 3, 4);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left 120 (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(-0.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 1 X3, 1 X2 } { 1 Y3 };\n", \ "\n", \ "newlines\n", \ "(3, 4);\n", \ "(4, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.866) { 2 X2, 1 X3 } { 1 Y3 };\n", \ "(1, 1.732) { 2 X2, 2 X3 } { 2 Y3 };\n", \ "(0, 1.732) { 1 X2, 2 X3 } { 2 Y3 };\n", \ "(-0.5, 0.866) { 1 X3 } {1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 4);\n", \ "(1, 4, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Left Right 120 (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(-0.5, 0.866);\n", \ "(1.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 1) del;\n", \ "(2, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n", \ "\n", \ "newlines\n", \ "(3, 5);\n", \ "(5, 4);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.5, 0.866) { 1 X4 } { 1 Y4 };\n", \ "(1, 1.299) { -0.5 X2, 0.375 X3, 1.125 X4 } { -0.5 Y2, 0.375 Y3, 1.125 Y4 };\n", \ "(0, 1.299) { 1.125 X3, 0.375 X4 } { 1.125 Y3, 0.375 Y4 };\n", \ "(-0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 5);\n", \ "(3, 1, 5);\n", \ "(2, 4, 5);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "rule \"Fill Triangle\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(2, 3) del;\n", \ "(3, 1) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { 1 Y2 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"Vis A Vis (1)\"\n", \ "\n", \ "quality 1\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(0.5, 0.866);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "\n", \ "newpoints\n", \ "\n", \ "newlines\n", \ "(1, 3);\n", \ "(3, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(1.2, 0.693) { 0.8 X2, 0.8 X3 } { 0.8 Y2, 0.8 Y3 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(-0.2, 0.693) { -0.6 X2, 0.8 X3 } { -0.6 Y2, 0.8 Y3 };\n", \ "\n", \ "freearea2\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { };\n", \ "(0.75, 0.433) { 0.5 X2, 0.5 X3 } { 0.5 Y2, 0.5 Y3 };\n", \ "(0.5, 0.866) { 1 X3 } { 1 Y3 };\n", \ "(0.25, 0.433) { 0.5 X3 } { 0.5 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 3);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "rule \"2 h Vis A Vis (1)\"\n", \ "\n", \ "quality 3\n", \ "\n", \ "mappoints\n", \ "(0, 0);\n", \ "(1, 0);\n", \ "(1, 1.732);\n", \ "(0, 1.732);\n", \ "\n", \ "maplines\n", \ "(1, 2) del;\n", \ "(3, 4) del;\n", \ "\n", \ "newpoints\n", \ "(0.5, 0.866) { 0.25 X2, 0.25 X3, 0.25 X4 } { 0.25 Y2, 0.25 Y3, 0.25 Y4 };\n", \ "\n", \ "newlines\n", \ "(1, 5);\n", \ "(5, 4);\n", \ "(3, 5);\n", \ "(5, 2);\n", \ "\n", \ "freearea\n", \ "(0, 0);\n", \ "(1, 0) { 1 X2 } { 1 Y2 };\n", \ "(1.5, 0.866) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y2, 0.75 Y3, -0.25 Y4 };\n", \ "(1, 1.732) { 1 X3 } { 1 Y3 };\n", \ "(0, 1.732) { 1 X4 } { 1 Y4 };\n", \ "(-0.5, 0.866) { 0.75 X4, -0.25 X2, -0.25 X3 } { 0.75 Y4, -0.25 Y3 };\n", \ "\n", \ "elements\n", \ "(1, 2, 5);\n", \ "(3, 4, 5);\n", \ "\n", \ "endrule\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ "\n", \ 0};

	  public static void GetPureBadness(Mesh mesh, ref Array<double> pure_badness, BitArray isnewpoint)
	  {
		//const int ne = mesh.GetNE();
		int np = mesh.GetNP();

		pure_badness.SetSize(np + PointIndex.BASE+1);
		pure_badness = -1;

		Array< Point < 3> > backup = new Array< Point < 3> >(np);

		for (int i = 0; i < np; i++)
		{
		backup[i] = new Point < 3>(new mesh.Point(i + 1));

		if (isnewpoint.Test(i + PointIndex.BASE) && mesh.mlbetweennodes[i + PointIndex.BASE][0] > 0)
		{
			mesh.Point(i + 1) = Center(new mesh.Point(mesh.mlbetweennodes[i + PointIndex.BASE][0]), new mesh.Point(mesh.mlbetweennodes[i + PointIndex.BASE][1]));
		}
		}
		for (ElementIndex i = 0; i < mesh.GetNE(); i++)
		{
		double bad = mesh[i].CalcJacobianBadness(mesh.Points());
		for (int j = 0; j < mesh[i].GetNP(); j++)
		{
		  if (bad > pure_badness[mesh[i][j]])
		  {
			pure_badness[mesh[i][j]] = bad;
		  }
		}

		// save maximum
		if (bad > pure_badness.Last())
		{
		  pure_badness.Last() = bad;
		}
		}

		for (int i = 0; i < np; i++)
		{
		mesh.Point(i + 1) = *backup[i];
		backup[i] = null;
		}
	  }

	  public static double Validate(Mesh mesh, Array<ElementIndex> bad_elements, Array<double> pure_badness, double max_worsening, bool uselocalworsening, Array<double> quality_loss = null)
	  {
		PrintMessage(3, "!!!! Validating !!!!");
		//if(max_worsening > 0)
		//  (*testout) << "badness " << counter++ << endl;

		bad_elements.SetSize(0);

		double loc_pure_badness = -1;

		if (!uselocalworsening)
		{
		  loc_pure_badness = pure_badness.Last(); // maximum is saved at last position
		}


		double worsening = -1;
		ElementIndex ind = new ElementIndex();

		if (quality_loss != null)
		{
		  quality_loss.SetSize(mesh.GetNE());
		}

		for (ElementIndex i = 0; i < mesh.GetNE(); i++)
		{
		if (uselocalworsening)
		{
			loc_pure_badness = -1;
			for (int j = 0; j < mesh[i].GetNP(); j++)
			{
			  if (pure_badness[mesh[i][j]] > loc_pure_badness)
			  {
			loc_pure_badness = pure_badness[mesh[i][j]];
			  }
			}
		}


		double bad = mesh[i].CalcJacobianBadness(mesh.Points());
		if (bad > 1e10 || (max_worsening > 0 && bad > loc_pure_badness * max_worsening))
		{
		  bad_elements.Append(i);
		}


		if (max_worsening > 0)
		{
			double actw = bad / loc_pure_badness;
			if (quality_loss != null)
			{
			  quality_loss[i] = actw;
			}

			if (actw > worsening)
			{
			worsening = actw;
			ind.CopyFrom(i);
			}
		}
		}
		return worsening;
	  }

	  public static void RepairBisection(Mesh mesh, Array<ElementIndex> bad_elements, BitArray isnewpoint, Refinement refinement, Array<double> pure_badness, double max_worsening, bool uselocalworsening, Array< Array<int,PointIndex.BASE> > idmaps)
	  {
		ostringstream ostrstr = new ostringstream();

		const int maxtrials = 100;

		//bool doit;
		//cout << "DOIT: " << flush;
		//cin >> doit;

		int ne = mesh.GetNE();
		int np = mesh.GetNP();

		int numbadneighbours = 3;
		const int numtopimprove = 3;

		PrintMessage(1, "repairing");

		PushStatus("Repair Bisection");

		Array<Point < 3> > should = new Array<Point < 3> >(np);
		Array<Point < 3> > can = new Array<Point < 3> >(np);
		Array<Vec < 3> > nv = new Array<Vec < 3> >(np);
		for (int i = 0; i < np; i++)
		{
		nv[i] = new Vec < 3>;
		should[i] = new Point < 3>;
		can[i] = new Point < 3>;
		}

		BitArray isboundarypoint = new BitArray(np);
		BitArray isedgepoint = new BitArray(np);
		isboundarypoint.Clear();
		isedgepoint.Clear();

		for (int i = 1; i <= mesh.GetNSeg(); i++)
		{
		Segment seg = mesh.LineSegment(i);
		isedgepoint.Set(new netgen.Segment(seg[0]));
		isedgepoint.Set(new netgen.Segment(seg[1]));
		}

		Array<int> surfaceindex = new Array<int>(np);
		surfaceindex = -1;

		for (int i = 1; i <= mesh.GetNSE(); i++)
		{
		Element2d sel = mesh.SurfaceElement(i);
		for (int j = 1; j <= sel.GetNP(); j++)
		{
		  if (!isedgepoint.Test(sel.PNum(j)))
		  {
			  isboundarypoint.Set(sel.PNum(j));
			  surfaceindex[sel.PNum(j) - PointIndex.BASE] = mesh.GetFaceDescriptor(sel.GetIndex()).SurfNr();
		  }
		}
		}



		Validate(mesh, bad_elements, pure_badness, ((uselocalworsening) ? (0.8 * (max_worsening - 1.0) + 1.0) : (0.1 * (max_worsening - 1.0) + 1.0)), uselocalworsening); // -> larger working area
		BitArray working_elements = new BitArray(ne);
		BitArray working_points = new BitArray(np);

		GetWorkingArea(working_elements, working_points, mesh, bad_elements, numbadneighbours);
		//working_elements.Set();
		//working_points.Set();

		ostrstr.str("");
		ostrstr << "worsening: " << Validate(mesh, bad_elements, pure_badness, max_worsening, uselocalworsening);
		PrintMessage(4, ostrstr.str());



		int auxnum = 0;
		for (int i = 1; i <= np; i++)
		{
		  if (working_points.Test(i))
		  {
		auxnum++;
		  }
		}

		ostrstr.str("");
		ostrstr << "Percentage working points: " << 100.0 * (double)auxnum / np;
		PrintMessage(5, ostrstr.str());


		BitArray isworkingboundary = new BitArray(np);
		for (int i = 1; i <= np; i++)
		{
		  if (working_points.Test(i) && isboundarypoint.Test(i))
		  {
		isworkingboundary.Set(i);
		  }
		  else
		  {
		isworkingboundary.Clear(i);
		  }
		}


		for (int i = 0; i < np; i++)
		{
		  *should[i] = new mesh.Point(i + 1);
		}


		for (int i = 0; i < np; i++)
		{
		if (isnewpoint.Test(i + PointIndex.BASE) && mesh.mlbetweennodes[i + PointIndex.BASE][0] > 0)
		{
		  *can[i] = Center(can[mesh.mlbetweennodes[i + PointIndex.BASE][0] - PointIndex.BASE], can[mesh.mlbetweennodes[i + PointIndex.BASE][1] - PointIndex.BASE]);
		}
		else
		{
		  *can[i] = new mesh.Point(i + 1);
		}
		}


		int cnttrials = 1;

		double lamedge = 0.5;
		double lamface = 0.5;

		double facokedge = 0;
		double facokface = 0;
		double factryedge;
		double factryface = 0;

		double oldlamedge;
		double oldlamface;

		MeshOptimize2d optimizer2d = refinement.Get2dOptimizer();
		if (optimizer2d == null)
		{
		cerr << "No 2D Optimizer!" << "\n";
		return;
		}

		while ((facokedge < 1.0 - 1e-8 || facokface < 1.0 - 1e-8) && cnttrials < maxtrials && multithread.terminate != 1)
		{
		(*testout) << "   facokedge " << facokedge << " facokface " << facokface << " cnttrials " << cnttrials << "\n" << " perc. " << 95.0 * max2(min2(facokedge, facokface), (double)cnttrials / (double)maxtrials) << "\n";

		SetThreadPercent(95.0 * max2(min2(facokedge, facokface), (double)cnttrials / (double)maxtrials));

		ostrstr.str("");
		ostrstr << "max. worsening " << max_worsening;
		PrintMessage(5, ostrstr.str());
		oldlamedge = lamedge;
		lamedge *= 6;
		if (lamedge > 2)
		{
		  lamedge = 2;
		}

		if (1 == 1 || facokedge < 1.0 - 1e-8)
		{
			for (int i = 0; i < nv.Size(); i++)
			{
			  *nv[i] = Vec < 3>(0,0,0);
			}
			for (int i = 1; i <= mesh.GetNSE(); i++)
			{
			Element2d sel = mesh.SurfaceElement(i);
			Vec < 3> auxvec = Cross(new mesh.Point(sel.PNum(2)) - new mesh.Point(sel.PNum(1)), new mesh.Point(sel.PNum(3)) - new mesh.Point(sel.PNum(1)));
			auxvec.Normalize();
			for (int j = 1; j <= sel.GetNP(); j++)
			{
			  if (!isedgepoint.Test(sel.PNum(j)))
			  {
				*nv[sel.PNum(j) - PointIndex.BASE] += auxvec;
			  }
			}
			}
			for (int i = 0; i < nv.Size(); i++)
			{
			  nv[i].Normalize();
			}


			do // move edges
			{
			lamedge *= 0.5;
			cnttrials++;
			if (cnttrials % 10 == 0)
			{
			  max_worsening *= 1.1;
			}


			factryedge = lamedge + (1.0 - lamedge) * facokedge;

			ostrstr.str("");
			ostrstr << "lamedge = " << lamedge << ", trying: " << factryedge;
			PrintMessage(5, ostrstr.str());


			for (int i = 1; i <= np; i++)
			{
				if (isedgepoint.Test(i))
				{
				for (int j = 0; j < 3; j++)
				{
				  mesh.Point(i)(j) = lamedge * should.Get(i)(j) + (1.0 - lamedge) * can.Get(i)(j);
				}
				}
				else
				{
				  mesh.Point(i) = *can.Get(i);
				}
			}
			if (facokedge < 1.0 - 1e-8)
			{
				ostrstr.str("");
				ostrstr << "worsening: " << Validate(mesh, bad_elements, pure_badness, max_worsening, uselocalworsening);

				PrintMessage(5, ostrstr.str());
			}
			else
			{
			  Validate(mesh, bad_elements, pure_badness, -1, uselocalworsening);
			}


			ostrstr.str("");
			ostrstr << bad_elements.Size() << " bad elements";
			PrintMessage(5, ostrstr.str());
			} while (bad_elements.Size() > 0 && cnttrials < maxtrials && multithread.terminate != 1);
		}

		if (cnttrials < maxtrials && multithread.terminate != 1)
		{
			facokedge = factryedge;

			// smooth faces
			mesh.CalcSurfacesOfNode();

			MeshingParameters dummymp = new MeshingParameters();
			mesh.ImproveMeshJacobianOnSurface(dummymp, isworkingboundary, nv, OPTIMIZEGOAL.OPT_QUALITY, idmaps);

			for (int i = 1; i <= np; i++)
			{
			  *can.Elem(i) = new mesh.Point(i);
			}

			if (optimizer2d != null)
			{
			  optimizer2d.ProjectBoundaryPoints(surfaceindex, can, should);
			}
		}


		oldlamface = lamface;
		lamface *= 6;
		if (lamface > 2)
		{
		  lamface = 2;
		}


		if (cnttrials < maxtrials && multithread.terminate != 1)
		{

			do // move faces
			{
			lamface *= 0.5;
			cnttrials++;
			if (cnttrials % 10 == 0)
			{
			  max_worsening *= 1.1;
			}
			factryface = lamface + (1.0 - lamface) * facokface;

			ostrstr.str("");
			ostrstr << "lamface = " << lamface << ", trying: " << factryface;
			PrintMessage(5, ostrstr.str());


			for (int i = 1; i <= np; i++)
			{
				if (isboundarypoint.Test(i))
				{
				for (int j = 0; j < 3; j++)
				{
				  mesh.Point(i)(j) = lamface * should.Get(i)(j) + (1.0 - lamface) * can.Get(i)(j);
				}
				}
				else
				{
				  mesh.Point(i) = *can.Get(i);
				}
			}

			ostrstr.str("");
			ostrstr << "worsening: " << Validate(mesh, bad_elements, pure_badness, max_worsening, uselocalworsening);
			PrintMessage(5, ostrstr.str());


			ostrstr.str("");
			ostrstr << bad_elements.Size() << " bad elements";
			PrintMessage(5, ostrstr.str());
			} while (bad_elements.Size() > 0 && cnttrials < maxtrials && multithread.terminate != 1);
		}



		if (cnttrials < maxtrials && multithread.terminate != 1)
		{
			facokface = factryface;
			// smooth interior

			mesh.CalcSurfacesOfNode();

			MeshingParameters dummymp = new MeshingParameters();
			mesh.ImproveMeshJacobian(dummymp, OPTIMIZEGOAL.OPT_QUALITY, working_points);
			//mesh.ImproveMeshJacobian (OPT_WORSTCASE,&working_points);


			for (int i = 1; i <= np; i++)
			{
			  *can.Elem(i) = new mesh.Point(i);
			}
		}

		//!
		if ((facokedge < 1.0 - 1e-8 || facokface < 1.0 - 1e-8) && cnttrials < maxtrials && multithread.terminate != 1)
		{
			MeshingParameters dummymp = new MeshingParameters();
			MeshOptimize3d optmesh = new MeshOptimize3d(dummymp);
			for (int i = 0; i < numtopimprove; i++)
			{
			optmesh.SwapImproveSurface(mesh, OPTIMIZEGOAL.OPT_QUALITY, working_elements, idmaps);
			optmesh.SwapImprove(mesh, OPTIMIZEGOAL.OPT_QUALITY, working_elements);

			}

			//	    mesh.mglevels = 1;


			ne = mesh.GetNE();
			working_elements.SetSize(ne);


			for (int i = 1; i <= np; i++)
			{
			  mesh.Point(i) = *should.Elem(i);
			}

			Validate(mesh, bad_elements, pure_badness, ((uselocalworsening) ? (0.8 * (max_worsening - 1.0) + 1.0) : (0.1 * (max_worsening - 1.0) + 1.0)), uselocalworsening);

			if (lamedge < oldlamedge || lamface < oldlamface)
			{
			  numbadneighbours++;
			}
			GetWorkingArea(working_elements, working_points, mesh, bad_elements, numbadneighbours);
			for (int i = 1; i <= np; i++)
			{
			  if (working_points.Test(i) && isboundarypoint.Test(i))
			  {
			isworkingboundary.Set(i);
			  }
			  else
			  {
			isworkingboundary.Clear(i);
			  }
			}
			auxnum = 0;
			for (int i = 1; i <= np; i++)
			{
			  if (working_points.Test(i))
			  {
			auxnum++;
			  }
			}


			ostrstr.str("");
			ostrstr << "Percentage working points: " << 100.0 * (double)auxnum / np;
			PrintMessage(5, ostrstr.str());

			for (int i = 1; i <= np; i++)
			{
			  mesh.Point(i) = *can.Elem(i);
			}
		}
		//!

		}

		MeshingParameters dummymp = new MeshingParameters();
		MeshOptimize3d optmesh = new MeshOptimize3d(dummymp);
		for (int i = 0; i < numtopimprove && multithread.terminate != 1; i++)
		{
		optmesh.SwapImproveSurface(mesh, OPTIMIZEGOAL.OPT_QUALITY, null, idmaps);
		optmesh.SwapImprove(mesh, OPTIMIZEGOAL.OPT_QUALITY);
		//mesh.UpdateTopology();
		}
		mesh.UpdateTopology();
		/*
		if(cnttrials < 100)
		  {
		nv = Vec3d(0,0,0);
		for (int i = 1; i <= mesh.GetNSE(); i++)
		  {
			const Element2d & sel = mesh.SurfaceElement(i);
			Vec3d auxvec = Cross(mesh.Point(sel.PNum(2))-mesh.Point(sel.PNum(1)),
					 mesh.Point(sel.PNum(3))-mesh.Point(sel.PNum(1)));
			auxvec.Normalize();
			for (int j = 1; j <= sel.GetNP(); j++)
			  if(!isedgepoint.Test(sel.PNum(j)))
			nv[sel.PNum(j) - PointIndex::BASE] += auxvec;
		  }
		for(int i=0; i<nv.Size(); i++)
		  nv[i].Normalize();
	
	
		mesh.ImproveMeshJacobianOnSurface(isboundarypoint,nv,OPT_QUALITY);
		mesh.CalcSurfacesOfNode();
			// smooth interior
	
	
		for (int i = 1; i <= np; i++)
		  if(isboundarypoint.Test(i))
			can.Elem(i) = mesh.Point(i);
	
		if(optimizer2d)
		  optimizer2d->ProjectBoundaryPoints(surfaceindex,can,should);
	
	
		for (int i = 1; i <= np; i++)
		  if(isboundarypoint.Test(i))
			for(int j=1; j<=3; j++)
			  mesh.Point(i).X(j) = should.Get(i).X(j);
		  }
		*/


		if (cnttrials == maxtrials)
		{
		for (int i = 1; i <= np; i++)
		{
		  mesh.Point(i) = *should.Get(i);
		}

		Validate(mesh, bad_elements, pure_badness, max_worsening, uselocalworsening);

		for (int i = 0; i < bad_elements.Size(); i++)
		{
			ostrstr.str("");
			ostrstr << "bad element:" << "\n" << mesh[bad_elements[i]][0] << ": " << new mesh.Point(mesh[bad_elements[i]][0]) << "\n" << mesh[bad_elements[i]][1] << ": " << new mesh.Point(mesh[bad_elements[i]][1]) << "\n" << mesh[bad_elements[i]][2] << ": " << new mesh.Point(mesh[bad_elements[i]][2]) << "\n" << mesh[bad_elements[i]][3] << ": " << new mesh.Point(mesh[bad_elements[i]][3]);
			PrintMessage(5, ostrstr.str());
		}
		for (int i = 1; i <= np; i++)
		{
		  mesh.Point(i) = *can.Get(i);
		}
		}

		for (int i = 0; i < np; i++)
		{
		nv[i] = null;
		can[i] = null;
		should[i] = null;
		}

		PopStatus();
	  }


	  public static void GetWorkingArea(BitArray working_elements, BitArray working_points, Mesh mesh, Array<ElementIndex> bad_elements, int width)
	  {
		working_elements.Clear();
		working_points.Clear();

		for (int i = 0; i < bad_elements.Size(); i++)
		{
		working_elements.Set(bad_elements[i]);
		Element el = mesh[bad_elements[i]];
		for (int j = 1; j <= el.GetNP(); j++)
		{
		  working_points.Set(el.PNum(j));
		}
		}


		for (int i = 0; i < width; i++)
		{
		for (ElementIndex j = 0; j < mesh.GetNE(); j++)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: if(!working_elements.Test(j))
			if (!working_elements.Test(new netgen.ElementIndex(j)))
			{
			Element el = mesh[j];
			bool set_active = false;

			for (int k = 1; !set_active && k <= el.GetNP(); k++)
			{
			  set_active = working_points.Test(el.PNum(k));
			}

			if (set_active)
			{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: working_elements.Set(j);
			  working_elements.Set(new netgen.ElementIndex(j));
			}
			}
		}

		for (ElementIndex j = 0; j < mesh.GetNE(); j++)
		{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: if(working_elements.Test(j))
			if (working_elements.Test(new netgen.ElementIndex(j)))
			{
			Element el = mesh[j];
			for (int k = 1; k <= el.GetNP(); k++)
			{
			  working_points.Set(el.PNum(k));
			}
			}
		}
		}
	  }
	}
}