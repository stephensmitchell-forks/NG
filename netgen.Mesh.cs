namespace netgen
{

	public class Mesh
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void BuildCurvedElements(Refinement @ref, int aorder, bool arational)
		  {
			GetCurvedElements().BuildCurvedElements(@ref, aorder, arational);
        
			for (SegmentIndex seg = 0; seg < GetNSeg(); seg++)
			{
			  this[seg].SetCurved(GetCurvedElements().IsSegmentCurved(seg));
			}
			for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			{
			  this[sei].SetCurved(GetCurvedElements().IsSurfaceElementCurved(sei));
			}
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			{
			  this[ei].SetCurved(GetCurvedElements().IsElementCurved(ei));
			}
        
			SetNextMajorTimeStamp();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void SendRecvMesh()
		  {
			int id = GetCommunicator().Rank();
			int np = GetCommunicator().Size();
        
			if (np == 1)
			{
			  throw new Exception("SendRecvMesh called, but only one rank in communicator!!");
			}
        
			if (id == 0)
			{
			  PrintMessage(1, "Send/Receive mesh");
			}
        
			// Why is this here??
			if (id == 0)
			{
			paralleltop.SetNV(GetNV());
			paralleltop.SetNE(GetNE());
			paralleltop.SetNSegm(GetNSeg());
			paralleltop.SetNSE(GetNSE());
			}
        
			if (id == 0)
			{
			  SendMesh();
			}
			else
			{
			  ReceiveParallelMesh();
			}
        
			paralleltop.UpdateCoarseGrid();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void SendMesh()
		  {
			Array<MPI_Request> sendrequests = new Array<MPI_Request>();
        
			NgMPI_Comm comm = GetCommunicator();
			int id = comm.Rank();
			int ntasks = comm.Size();
        
			int dim = GetDimension();
			MyMPI_Bcast(dim, new NgMPI_Comm(comm));
        
        
			// If the topology is not already updated, we do not need to
			// build edges/faces.
		//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
			var top = const_cast<MeshTopology&>(GetTopology());
			if (top.NeedsUpdate())
			{
			  top.SetBuildEdges(false);
			  top.SetBuildFaces(false);
			  top.Update();
			}
        
			PrintMessage(3, "Sending nr of elements");
        
			Array<int> num_els_on_proc = new Array<int>(ntasks);
			num_els_on_proc = 0;
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			{
			  // num_els_on_proc[(*this)[ei].GetPartition()]++;
			  num_els_on_proc[vol_partition[ei]]++;
			}
        
			MPI_Scatter(num_els_on_proc[0], 1, MPI_INT, MPI_IN_PLACE, -1, MPI_INT, 0, comm);
        
			TABLE<ElementIndex> els_of_proc = new TABLE<ElementIndex>(num_els_on_proc);
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			{
			  // els_of_proc.Add ( (*this)[ei].GetPartition(), ei);
			  els_of_proc.Add(vol_partition[ei], ei);
			}
        
			PrintMessage(3, "Building vertex/proc mapping");
        
			Array<int> num_sels_on_proc = new Array<int>(ntasks);
			num_sels_on_proc = 0;
			for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
			{
			  // num_sels_on_proc[(*this)[ei].GetPartition()]++;
			  num_sels_on_proc[surf_partition[ei]]++;
			}
        
			TABLE<SurfaceElementIndex> sels_of_proc = new TABLE<SurfaceElementIndex>(num_sels_on_proc);
			for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
			{
			  // sels_of_proc.Add ( (*this)[ei].GetPartition(), ei);
			  sels_of_proc.Add(surf_partition[ei], ei);
			}
        
			Array<int> num_segs_on_proc = new Array<int>(ntasks);
			num_segs_on_proc = 0;
			for (SegmentIndex ei = 0; ei < GetNSeg(); ei++)
			{
			  // num_segs_on_proc[(*this)[ei].GetPartition()]++;
			  num_segs_on_proc[seg_partition[ei]]++;
			}
        
			TABLE<SegmentIndex> segs_of_proc = new TABLE<SegmentIndex>(num_segs_on_proc);
			for (SegmentIndex ei = 0; ei < GetNSeg(); ei++)
			{
			  segs_of_proc.Add(seg_partition[ei], ei);
			}
        
        
			/**
			        ----- STRATEGY FOR PERIODIC MESHES -----
		
			   Whenever two vertices are identified by periodicity, any proc 
			   that gets one of the vertices actually gets both of them.
			   This has to be transitive, that is, if
			   a <-> b and  b <-> c,
			   then any proc that has vertex a also has vertices b and c!
		
			   Surfaceelements and Segments that are identified by
			   periodicity are treated the same way.
			   
			   We need to duplicate these so we have containers to
			   hold the edges/facets. Afaik, a mesh cannot have nodes 
			   that are not part of some sort of element.
		
			 **/
        
			/** First, we build tables for vertex identification. **/
			Array<INDEX_2> per_pairs = new Array<INDEX_2>();
			Array<INDEX_2> pp2 = new Array<INDEX_2>();
			var idents = GetIdentifications();
			bool has_periodic = false;
			for (int idnr = 1; idnr < idents.GetMaxNr() + 1; idnr++)
			{
			if (idents.GetType(idnr) != Identifications.PERIODIC)
			{
				continue;
			}
			has_periodic = true;
			idents.GetPairs(idnr, pp2);
			per_pairs.Append(pp2);
			}
			Array<int, PointIndex.BASE> npvs = new Array<int, PointIndex.BASE>(GetNV());
			npvs = 0;
			for (int k = 0; k < per_pairs.Size(); k++)
			{
			  npvs[per_pairs[k].I1()]++;
			  npvs[per_pairs[k].I2()]++;
			}
        
			/** for each vertex, gives us all identified vertices **/
			TABLE<PointIndex, PointIndex.BASE> per_verts = new TABLE<PointIndex, PointIndex.BASE>(npvs);
			for (int k = 0; k < per_pairs.Size(); k++)
			{
			  per_verts.Add(per_pairs[k].I1(), per_pairs[k].I2());
			  per_verts.Add(per_pairs[k].I2(), per_pairs[k].I1());
			}
			for (int k = PointIndex.BASE; k < GetNV() + PointIndex.BASE; k++)
			{
			  BubbleSort(per_verts[k]);
			}
        
			/** The same table as per_verts, but TRANSITIVE!! **/
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_per_verts_trans = [&](auto f)
			var iterate_per_verts_trans = (f UnnamedParameter) =>
			{
			  Array<int> allvs = new Array<int>();
			  for (int k = PointIndex.BASE; k < GetNV() + PointIndex.BASE; k++)
			  {
			  allvs.SetSize(0);
			  allvs.Append(per_verts[k]);
			  bool changed = true;
			  while (changed)
			  {
				changed = false;
				for (int j = 0; j < allvs.Size(); j++)
				{
				var pervs2 = per_verts[allvs[j]];
				for (int l = 0; l < pervs2.Size(); l++)
				{
					var addv = pervs2[l];
					if (allvs.Contains(addv) || addv == k)
					{
						continue;
					}
					changed = true;
					allvs.Append(addv);
				}
				}
			  }
			  f(k, allvs);
			  }
			};
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_per_verts_trans([&](auto k, auto & allvs)
			iterate_per_verts_trans((k UnnamedParameter, allvs) =>
			{
			npvs[k] = allvs.Size();
			});
			TABLE<PointIndex, PointIndex.BASE> per_verts_trans = new TABLE<PointIndex, PointIndex.BASE>(npvs);
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_per_verts_trans([&](auto k, auto & allvs)
			iterate_per_verts_trans((k UnnamedParameter, allvs) =>
			{
			for (int j = 0; j < allvs.Size(); j++)
			{
			  per_verts_trans.Add(k, allvs[j]);
			}
			});
			for (int k = PointIndex.BASE; k < GetNV() + PointIndex.BASE; k++)
			{
			  BubbleSort(per_verts_trans[k]);
			}
        
			/** Now we build the vertex-data to send to the workers. **/
			Array<int, PointIndex.BASE> vert_flag = new Array<int, PointIndex.BASE>(GetNV());
			Array<int, PointIndex.BASE> num_procs_on_vert = new Array<int, PointIndex.BASE>(GetNV());
			Array<int> num_verts_on_proc = new Array<int>(ntasks);
			num_verts_on_proc = 0;
			num_procs_on_vert = 0;
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_vertices = [&](auto f)
			var iterate_vertices = (f UnnamedParameter) =>
			{
			  vert_flag = -1;
			  for (int dest = 1; dest < ntasks; dest++)
			  {
			  FlatArray<ElementIndex> els = els_of_proc[dest];
			  for (int hi = 0; hi < els.Size(); hi++)
			  {
				  Element el = this[els[hi]];
				  for (int i = 0; i < el.GetNP(); i++)
				  {
				f(el[i], dest);
				  }
			  }
			  FlatArray<SurfaceElementIndex> sels = sels_of_proc[dest];
			  for (int hi = 0; hi < sels.Size(); hi++)
			  {
				  Element2d el = this[sels[hi]];
				  for (int i = 0; i < el.GetNP(); i++)
				  {
				f(el[i], dest);
				  }
			  }
			  FlatArray<SegmentIndex> segs = segs_of_proc[dest];
			  for (int hi = 0; hi < segs.Size(); hi++)
			  {
				  Segment el = this[segs[hi]];
				  for (int i = 0; i < 2; i++)
				  {
				f(el[i], dest);
				  }
			  }
			  }
			};
			/** count vertices per proc and procs per vertex **/
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_vertices([&](auto vertex, auto dest)
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
			iterate_vertices((vertex UnnamedParameter, dest UnnamedParameter2) =>
			{
			var countit = (vertex UnnamedParameter, dest UnnamedParameter2) =>
			{
			  if (vert_flag[vertex] < dest)
			  {
				  vert_flag[vertex] = dest;
				  num_verts_on_proc[dest]++;
				  num_procs_on_vert[vertex]++;
				  GetParallelTopology().SetDistantPNum(dest, vertex);
			  }
			};
			countit(vertex, dest);
			var pers = per_verts_trans[vertex];
			for (int j = 0; j < pers.Size(); j++)
			{
			  countit(pers[j], dest);
			}
			});
			TABLE<PointIndex> verts_of_proc = new TABLE<PointIndex>(num_verts_on_proc);
			TABLE<int, PointIndex.BASE> procs_of_vert = new TABLE<int, PointIndex.BASE>(num_procs_on_vert);
			TABLE<int, PointIndex.BASE> loc_num_of_vert = new TABLE<int, PointIndex.BASE>(num_procs_on_vert);
			/** Write vertex/proc mappingfs to tables **/
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_vertices([&](auto vertex, auto dest)
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
			iterate_vertices((vertex UnnamedParameter, dest UnnamedParameter2) =>
			{
			var addit = (vertex UnnamedParameter, dest UnnamedParameter2) =>
			{
			  if (vert_flag[vertex] < dest)
			  {
				  vert_flag[vertex] = dest;
				  procs_of_vert.Add(vertex, dest);
			  }
			};
			addit(vertex, dest);
			var pers = per_verts_trans[vertex];
			for (int j = 0; j < pers.Size(); j++)
			{
			  addit(pers[j], dest);
			}
			});
			/** 
			local vertex numbers on distant procs 
			(I think this was only used for debugging??) 
			**/
			for (int vert = 1; vert <= GetNP(); vert++)
			{
			FlatArray<int> procs = procs_of_vert[vert];
			for (int j = 0; j < procs.Size(); j++)
			{
				int dest = procs[j];
				verts_of_proc.Add(dest, vert);
				loc_num_of_vert.Add(vert, verts_of_proc[dest].Size());
			}
			}
			PrintMessage(3, "Sending Vertices - vertices");
        
			for (int dest = 1; dest < ntasks; dest++)
			{
			FlatArray<PointIndex> verts = verts_of_proc[dest];
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: sendrequests.Append(MyMPI_ISend(verts, dest, MPI_TAG_MESH+1, comm));
			sendrequests.Append(MyMPI_ISend(new FlatArray<PointIndex>(verts), dest, MPI_TAG_MESH + 1, new NgMPI_Comm(comm)));
        
			int mptype = MeshPoint.MyGetMPIType();
        
			int numv = (int)verts.Size();
        
			int newtype;
			Array<int> blocklen = new Array<int>(numv);
			blocklen = 1;
        
		//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'reinterpret_cast' in C#:
			MPI_Type_indexed(numv, blocklen[0], reinterpret_cast<int> (verts[0]), mptype, newtype);
			MPI_Type_commit(newtype);
        
			MPI_Request request = new MPI_Request();
			MPI_Isend(points[0], 1, newtype, dest, MPI_TAG_MESH + 1, comm, request);
			sendrequests.Append(request);
			}
        
			Array<int> num_distpnums = new Array<int>(ntasks);
			num_distpnums = 0;
        
        
			/**
			   Next, we send the identifications themselfs.
			   
			   Info about periodic identifications sent to each proc is an array of
			   integers.
			   - maxidentnr
			   - type for each identification
			   - nr of pairs for each identification (each pair is local!)
			   - pairs for each periodic ident (global numbers)
			**/
			PrintMessage(3, "Sending Vertices - identifications");
			int maxidentnr = idents.GetMaxNr();
			Array<int> ppd_sizes = new Array<int>(ntasks);
			ppd_sizes = 1 + 2 * maxidentnr;
			for (int idnr = 1; idnr < idents.GetMaxNr() + 1; idnr++)
			{
			if (idents.GetType(idnr) != Identifications.PERIODIC)
			{
				continue;
			}
			idents.GetPairs(idnr, pp2);
			for (int j = 0; j < pp2.Size(); j++)
			{
				INDEX_2 pair = pp2[j];
				// both are on same procs!
				var ps = procs_of_vert[pair.I1()];
				for (int l = 0; l < ps.Size(); l++)
				{
				ppd_sizes[ps[l]] += 2;
				}
			}
			}
			TABLE<int> pp_data = new TABLE<int>(ppd_sizes);
			for (int dest = 0; dest < ntasks; dest++)
			{
			  pp_data.Add(dest, maxidentnr);
			}
			for (int dest = 0; dest < ntasks; dest++)
			{
			for (int idnr = 1; idnr < idents.GetMaxNr() + 1; idnr++)
			{
			  pp_data.Add(dest, idents.GetType(idnr));
			}
			for (int idnr = 1; idnr < idents.GetMaxNr() + 1; idnr++)
			{
			  pp_data.Add(dest, 0);
			}
			}
			for (int idnr = 1; idnr < idents.GetMaxNr() + 1; idnr++)
			{
			if (idents.GetType(idnr) != Identifications.PERIODIC)
			{
				continue;
			}
			idents.GetPairs(idnr, pp2);
			for (int j = 0; j < pp2.Size(); j++)
			{
				INDEX_2 pair = pp2[j];
				var ps = procs_of_vert[pair.I1()];
				for (int l = 0; l < ps.Size(); l++)
				{
				var p = ps[l];
				pp_data[p][maxidentnr + idnr]++;
				pp_data.Add(p, pair.I1());
				pp_data.Add(p, pair.I2());
				}
			}
			}
			Array<MPI_Request> req_per = new Array<MPI_Request>();
			for (int dest = 1; dest < ntasks; dest++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: req_per.Append(MyMPI_ISend(pp_data[dest], dest, MPI_TAG_MESH+1, comm));
			  req_per.Append(MyMPI_ISend(new TABLE<int>(pp_data[dest]), dest, MPI_TAG_MESH + 1, new NgMPI_Comm(comm)));
			}
			MPI_Waitall(req_per.Size(), req_per[0], MPI_STATUS_IGNORE);
        
			PrintMessage(3, "Sending Vertices - distprocs");
        
			for (int vert = 1; vert <= GetNP(); vert++)
			{
			FlatArray<int> procs = procs_of_vert[vert];
			for (int j = 0; j < procs.Size(); j++)
			{
			  num_distpnums[procs[j]] += 3 * (procs.Size() - 1);
			}
			}
        
			TABLE<int> distpnums = new TABLE<int>(num_distpnums);
        
			for (int vert = 1; vert <= GetNP(); vert++)
			{
			FlatArray<int> procs = procs_of_vert[vert];
			for (int j = 0; j < procs.Size(); j++)
			{
			  for (int k = 0; k < procs.Size(); k++)
			  {
				if (j != k)
				{
				distpnums.Add(procs[j], loc_num_of_vert[vert][j]);
				distpnums.Add(procs[j], procs_of_vert[vert][k]);
				distpnums.Add(procs[j], loc_num_of_vert[vert][k]);
				}
			  }
			}
			}
        
			for (int dest = 1; dest < ntasks; dest++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: sendrequests.Append(MyMPI_ISend(distpnums[dest], dest, MPI_TAG_MESH+1, comm));
			  sendrequests.Append(MyMPI_ISend(new TABLE<int>(distpnums[dest]), dest, MPI_TAG_MESH + 1, new NgMPI_Comm(comm)));
			}
        
        
        
			PrintMessage(3, "Sending elements");
        
			Array<int> elarraysize = new Array<int>(ntasks);
			elarraysize = 0;
			for (int ei = 1; ei <= GetNE(); ei++)
			{
			Element el = VolumeElement(ei);
			// int dest = el.GetPartition();
				int dest = vol_partition[ei - 1];
			elarraysize[dest] += 3 + el.GetNP();
			}
        
			TABLE<int> elementarrays = new TABLE<int>(elarraysize);
        
			for (int ei = 1; ei <= GetNE(); ei++)
			{
			Element el = VolumeElement(ei);
			// int dest = el.GetPartition();
				int dest = vol_partition[ei - 1];
        
			elementarrays.Add(dest, ei);
			elementarrays.Add(dest, el.GetIndex());
			elementarrays.Add(dest, el.GetNP());
			for (int i = 0; i < el.GetNP(); i++)
			{
			  elementarrays.Add(dest, el[i]);
			}
			}
        
			for (int dest = 1; dest < ntasks; dest++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: sendrequests.Append(MyMPI_ISend(elementarrays[dest], dest, MPI_TAG_MESH+2, comm));
			  sendrequests.Append(MyMPI_ISend(new TABLE<int>(elementarrays[dest]), dest, MPI_TAG_MESH + 2, new NgMPI_Comm(comm)));
			}
        
        
			PrintMessage(3, "Sending Face Descriptors");
        
			Array<double> fddata = new Array<double>(6 * GetNFD());
			for (int fdi = 1; fdi <= GetNFD(); fdi++)
			{
			fddata[6 * fdi - 6] = GetFaceDescriptor(fdi).SurfNr();
			fddata[6 * fdi - 5] = GetFaceDescriptor(fdi).DomainIn();
			fddata[6 * fdi - 4] = GetFaceDescriptor(fdi).DomainOut();
			fddata[6 * fdi - 3] = GetFaceDescriptor(fdi).BCProperty();
			fddata[6 * fdi - 2] = GetFaceDescriptor(fdi).domin_singular;
			fddata[6 * fdi - 1] = GetFaceDescriptor(fdi).domout_singular;
        
			}
			for (int dest = 1; dest < ntasks; dest++)
			{
			  sendrequests.Append(MyMPI_ISend(new Array<double>(fddata), dest, MPI_TAG_MESH + 3, new NgMPI_Comm(comm)));
			}
        
			/** Surface Elements **/
        
			PrintMessage(3, "Sending Surface elements");
			// build sel-identification
			uint nse = GetNSE();
			Array<SurfaceElementIndex> ided_sel = new Array<SurfaceElementIndex>(nse);
			ided_sel = -1;
			bool has_ided_sels = false;
			if (GetNE() && has_periodic) //we can only have identified surf-els if we have vol-els (right?)
			{
			Array<SurfaceElementIndex, 0> os1 = new Array<SurfaceElementIndex, 0>();
			Array<SurfaceElementIndex, 0> os2 = new Array<SurfaceElementIndex, 0>();
			for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			{
				if (ided_sel[sei] != -1)
				{
					continue;
				}
				Element2d sel = this[sei];
				FlatArray<PointIndex> points = sel.PNums();
				var ided1 = per_verts[points[0]];
				os1.SetSize(0);
				for (int j = 0; j < ided1.Size(); j++)
				{
				  os1.Append(GetTopology().GetVertexSurfaceElements(ided1[j]));
				}
				for (int j = 1; j < points.Size(); j++)
				{
				os2.SetSize(0);
				var p2 = points[j];
				var ided2 = per_verts[p2];
				for (int l = 0; l < ided2.Size(); l++)
				{
				  os2.Append(GetTopology().GetVertexSurfaceElements(ided2[l]));
				}
				for (int m = 0; m < os1.Size(); m++)
				{
				  if (!os2.Contains(os1[m]))
				  {
					os1.Delete(m);
					m--;
				  }
				}
				}
				if (!os1.Size())
				{
					continue;
				}
				if (os1.Size() > 1)
				{
				  throw new Exception("SurfaceElement identified with more than one other??");
				}
				Element2d sel2 = this[sei];
				FlatArray<PointIndex> points2 = sel2.PNums();
				has_ided_sels = true;
				ided_sel[sei] = os1[0];
				ided_sel[os1[0]] = sei;
			}
			}
			// build sel data to send
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_sels = [&](auto f)
			var iterate_sels = (f UnnamedParameter) =>
			{
			  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			  {
			  Element2d sel = this[sei];
			  // int dest = (*this)[sei].GetPartition();
				  int dest = surf_partition[sei];
			  f(sei, sel, dest);
			  if (ided_sel[sei] != -1)
			  {
				  // int dest2 = (*this)[ided_sel[sei]].GetPartition();
					  int dest2 = surf_partition[ided_sel[sei]];
				  f(sei, sel, dest2);
			  }
			  }
			};
			Array<int> nlocsel = new Array<int>(ntasks);
			Array<int> bufsize = new Array<int>(ntasks);
			nlocsel = 0;
			bufsize = 1;
			iterate_sels((SurfaceElementIndex sei, Element2d sel, int dest) =>
			{
			nlocsel[dest]++;
			bufsize[dest] += 4 + 2 * sel.GetNP();
			});
			TABLE<int> selbuf = new TABLE<int>(bufsize);
			for (int dest = 1; dest < ntasks; dest++)
			{
			  selbuf.Add(dest, nlocsel[dest]);
			}
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_sels([&](SurfaceElementIndex sei, const auto & sel, int dest)
			iterate_sels((SurfaceElementIndex sei, sel, int dest) =>
			{
			selbuf.Add(dest, sei);
			selbuf.Add(dest, sel.GetIndex());
			// selbuf.Add (dest, 0);
			selbuf.Add(dest, sel.GetNP());
			for (int ii = 1; ii <= sel.GetNP(); ii++)
			{
				selbuf.Add(dest, sel.PNum(ii));
				selbuf.Add(dest, sel.GeomInfoPi(ii).trignum);
			}
			});
			// distribute sel data
			for (int dest = 1; dest < ntasks; dest++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: sendrequests.Append(MyMPI_ISend(selbuf[dest], dest, MPI_TAG_MESH+4, comm));
			  sendrequests.Append(MyMPI_ISend(new TABLE<int>(selbuf[dest]), dest, MPI_TAG_MESH + 4, new NgMPI_Comm(comm)));
			}
        
        
			/** Segments **/
			PrintMessage(3, "Sending Edge Segments");
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_segs1 = [&](auto f)
			var iterate_segs1 = (f UnnamedParameter) =>
			{
			  Array<SegmentIndex> osegs1 = new Array<SegmentIndex>();
			  Array<SegmentIndex> osegs2 = new Array<SegmentIndex>();
			  Array<SegmentIndex> osegs_both = new Array<SegmentIndex>();
			  Array<int> type1 = new Array<int>();
			  Array<int> type2 = new Array<int>();
			  for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
			  {
			  Segment seg = this[segi];
			  int segnp = seg.GetNP();
			  PointIndex pi1 = seg[0];
			  var ided1 = per_verts[pi1];
			  PointIndex pi2 = seg[1];
			  var ided2 = per_verts[pi2];
			  if (!(ided1.Size() != 0 && ided2.Size() != 0))
			  {
				  continue;
			  }
			  osegs1.SetSize(0);
			  type1.SetSize(0);
			  for (int l = 0; l < ided1.Size(); l++)
			  {
				  var ospart = GetTopology().GetVertexSegments(ided1[l]);
				  for (int j = 0; j < ospart.Size(); j++)
				  {
				  if (osegs1.Contains(ospart[j]))
				  {
					throw new Exception("Periodic Mesh did something weird.");
				  }
				  osegs1.Append(ospart[j]);
				  type1.Append(idents.GetSymmetric(pi1, ided1[l]));
				  }
			  }
			  osegs2.SetSize(0);
			  type2.SetSize(0);
			  for (int l = 0; l < ided2.Size(); l++)
			  {
				  var ospart = GetTopology().GetVertexSegments(ided2[l]);
				  for (int j = 0; j < ospart.Size(); j++)
				  {
				  if (osegs2.Contains(ospart[j]))
				  {
					throw new Exception("Periodic Mesh did something weird.");
				  }
				  osegs2.Append(ospart[j]);
				  type2.Append(idents.GetSymmetric(pi2, ided2[l]));
				  }
			  }
			  osegs_both.SetSize(0);
			  for (int l = 0; l < osegs1.Size(); l++)
			  {
				var pos = osegs2.Pos(osegs1[l]);
				if (pos == -1)
				{
					continue;
				}
				if (type1[l] != type2[pos])
				{
					continue;
				}
				osegs_both.Append(osegs1[l]);
			  }
			  for (int l = 0; l < osegs_both.Size(); l++)
			  {
				int segnp2 = this[osegs_both[l]].GetNP();
				if (segnp != segnp2)
				{
				  throw new Exception("Tried to identify non-curved and curved Segment!");
				}
			  }
			  for (int l = 0; l < osegs_both.Size(); l++)
			  {
				f(segi, osegs_both[l]);
			  }
			  }
			};
			Array<int> per_seg_size = new Array<int>(GetNSeg());
			per_seg_size = 0;
			iterate_segs1((SegmentIndex segi1, SegmentIndex segi2) =>
			{
					  per_seg_size[segi1]++;
			});
			TABLE<SegmentIndex> per_seg = new TABLE<SegmentIndex>(per_seg_size);
			iterate_segs1((SegmentIndex segi1, SegmentIndex segi2) =>
			{
					  per_seg.Add(segi1, segi2);
			});
			// make per_seg transitive
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_per_seg_trans = [&](auto f)
			var iterate_per_seg_trans = (f UnnamedParameter) =>
			{
			  Array<SegmentIndex> allsegs = new Array<SegmentIndex>();
			  for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
			  {
			  allsegs.SetSize(0);
			  allsegs.Append(per_seg[segi]);
			  bool changed = true;
			  while (changed)
			  {
				  changed = false;
				  for (int j = 0; j < allsegs.Size(); j++)
				  {
				  var persegs2 = per_seg[allsegs[j]];
				  for (int l = 0; l < persegs2.Size(); l++)
				  {
					  var addseg = persegs2[l];
					  if (allsegs.Contains(addseg) || addseg == segi)
					  {
						  continue;
					  }
					  allsegs.Append(addseg);
					  changed = true;
				  }
				  }
			  }
			  f(segi, allsegs);
			  }
			};
			iterate_per_seg_trans((SegmentIndex segi, Array<SegmentIndex> segs) =>
			{
			for (int j = 0; j < segs.Size(); j++)
			{
			  per_seg_size[segi] = segs.Size();
			}
			});
			TABLE<SegmentIndex> per_seg_trans = new TABLE<SegmentIndex>(per_seg_size);
			iterate_per_seg_trans((SegmentIndex segi, Array<SegmentIndex> segs) =>
			{
			for (int j = 0; j < segs.Size(); j++)
			{
			  per_seg_trans.Add(segi, segs[j]);
			}
			});
			// build segment data
			Array<int> dests = new Array<int>();
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_segs2 = [&](auto f)
			var iterate_segs2 = (f UnnamedParameter) =>
			{
			for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
			{
				Segment seg = this[segi];
				dests.SetSize(0);
				// dests.Append(seg.GetPartition());
					dests.Append(seg_partition[segi]);
				for (int l = 0; l < per_seg_trans[segi].Size(); l++)
				{
				// int dest2 = (*this)[per_seg_trans[segi][l]].GetPartition();
						int dest2 = seg_partition[per_seg_trans[segi][l]];
				if (!dests.Contains(dest2))
				{
				  dests.Append(dest2);
				}
				}
				for (int l = 0; l < dests.Size(); l++)
				{
				  f(segi, seg, dests[l]);
				}
			}
			};
			Array<int> nloc_seg = new Array<int>(ntasks);
			// bufsize = 1; //was originally this - why??
			bufsize = 0;
			nloc_seg = 0;
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_segs2([&](auto segi, const auto & seg, int dest)
			iterate_segs2((segi UnnamedParameter, seg, int dest) =>
			{
					nloc_seg[dest]++;
					bufsize[dest] += 14;
			});
			TABLE<double> segm_buf = new TABLE<double>(bufsize);
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_segs2([&](auto segi, const auto & seg, int dest)
			iterate_segs2((segi UnnamedParameter, seg, int dest) =>
			{
					segm_buf.Add(dest, segi);
					segm_buf.Add(dest, seg.si);
					segm_buf.Add(dest, seg.pnums[0]);
					segm_buf.Add(dest, seg.pnums[1]);
					segm_buf.Add(dest, seg.geominfo[0].trignum);
					segm_buf.Add(dest, seg.geominfo[1].trignum);
					segm_buf.Add(dest, seg.surfnr1);
					segm_buf.Add(dest, seg.surfnr2);
					segm_buf.Add(dest, seg.edgenr);
					segm_buf.Add(dest, seg.epgeominfo[0].dist);
					segm_buf.Add(dest, seg.epgeominfo[1].edgenr);
					segm_buf.Add(dest, seg.epgeominfo[1].dist);
					segm_buf.Add(dest, seg.singedge_right);
					segm_buf.Add(dest, seg.singedge_left);
			});
			// distrubute segment data
			for (int dest = 1; dest < ntasks; dest++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: sendrequests.Append(MyMPI_ISend(segm_buf[dest], dest, MPI_TAG_MESH+5, comm));
			  sendrequests.Append(MyMPI_ISend(new TABLE<double>(segm_buf[dest]), dest, MPI_TAG_MESH + 5, new NgMPI_Comm(comm)));
			}
        
			PrintMessage(3, "now wait ...");
        
			MPI_Waitall(sendrequests.Size(), sendrequests[0], MPI_STATUS_IGNORE);
        
			PrintMessage(3, "Sending names");
        
			sendrequests.SetSize(3 * ntasks);
			/** Send bc/mat/cd*-names **/
			// nr of names
			int[] nnames = {0, 0, 0, 0};
			nnames[0] = materials.Size();
			nnames[1] = bcnames.Size();
			nnames[2] = GetNCD2Names();
			nnames[3] = GetNCD3Names();
			int tot_nn = nnames[0] + nnames[1] + nnames[2] + nnames[3];
			for (int k = 1; k < ntasks; k++)
			{
			   MPI_Isend(nnames, 4, MPI_INT, k, MPI_TAG_MESH + 6, comm, sendrequests[k]);
			}
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto iterate_names = [&](auto func)
			var iterate_names = (func UnnamedParameter) =>
			{
			  for (int k = 0; k < nnames[0]; k++)
			  {
				  func(materials[k]);
			  }
			  for (int k = 0; k < nnames[1]; k++)
			  {
				  func(bcnames[k]);
			  }
			  for (int k = 0; k < nnames[2]; k++)
			  {
				  func(cd2names[k]);
			  }
			  for (int k = 0; k < nnames[3]; k++)
			  {
				  func(cd3names[k]);
			  }
			};
			// sizes of names
			Array<int> name_sizes = new Array<int>(tot_nn);
			tot_nn = 0;
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_names([&](auto ptr)
			iterate_names((ptr UnnamedParameter) =>
			{
				name_sizes[tot_nn++] = (ptr == null) ? 0 : ptr.size();
			});
			for (int k = 1; k < ntasks; k++)
			{
			   MPI_Isend(name_sizes[0], tot_nn, MPI_INT, k, MPI_TAG_MESH + 6, comm, sendrequests[ntasks + k]);
			}
			// names
			int strs = 0;
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_names([&](auto ptr)
			iterate_names((ptr UnnamedParameter) =>
			{
				strs += (ptr == null) ? 0 : ptr.size();
			});
			Array<char> compiled_names = new Array<char>(strs);
			strs = 0;
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: iterate_names([&](auto ptr)
			iterate_names((ptr UnnamedParameter) =>
			{
			if (ptr == null)
			{
				return;
			}
			auto name = *ptr;
			for (int j = 0; j < name.size(); j++)
			{
				compiled_names[strs++] = name[j];
			}
			});
			for (int k = 1; k < ntasks; k++)
			{
			   MPI_Isend((compiled_names[0]), strs, MPI_CHAR, k, MPI_TAG_MESH + 6, comm, sendrequests[2 * ntasks + k]);
			}
        
			PrintMessage(3, "wait for names");
        
			MPI_Waitall(sendrequests.Size(), sendrequests[0], MPI_STATUS_IGNORE);
        
			MPI_Barrier(comm);
        
			PrintMessage(3, "send mesh complete");
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ReceiveParallelMesh()
		  {
			int timer = NgProfiler.CreateTimer("ReceiveParallelMesh");
			int timer_pts = NgProfiler.CreateTimer("Receive points");
			int timer_els = NgProfiler.CreateTimer("Receive elements");
			int timer_sels = NgProfiler.CreateTimer("Receive surface elements");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer);
        
			NgMPI_Comm comm = GetCommunicator();
			int id = comm.Rank();
			int ntasks = comm.Size();
        
			int dim;
			MyMPI_Bcast(dim, new NgMPI_Comm(comm));
			SetDimension(dim);
        
			// Receive number of local elements
			int nelloc;
			MPI_Scatter(null, 0, MPI_INT, nelloc, 1, MPI_INT, 0, comm);
			paralleltop.SetNE(nelloc);
        
			// string st;
        
			// receive vertices
			NgProfiler.StartTimer(timer_pts);
        
			Array<int> verts = new Array<int>();
			MyMPI_Recv(ref verts, 0, MPI_TAG_MESH + 1, new NgMPI_Comm(comm));
        
			int numvert = verts.Size();
			paralleltop.SetNV(numvert);
        
			// INDEX_CLOSED_HASHTABLE<int> glob2loc_vert_ht (3*numvert+1);
			INDEX_HASHTABLE<int> glob2loc_vert_ht = new INDEX_HASHTABLE<int>(3 * numvert + 1);
        
			for (int vert = 0; vert < numvert; vert++)
			{
			int globvert = verts[vert];
			paralleltop.SetLoc2Glob_Vert(vert + 1, globvert);
			glob2loc_vert_ht.Set(globvert, vert + 1);
			}
        
			for (int i = 0; i < numvert; i++)
			{
			  AddPoint(netgen.Point < 3> (0,0,0));
			}
        
			int mptype = MeshPoint.MyGetMPIType();
			MPI_Status status = new MPI_Status();
			MPI_Recv(points[1], numvert, mptype, 0, MPI_TAG_MESH + 1, comm, status);
        
			Array<int> pp_data = new Array<int>();
			MyMPI_Recv(ref pp_data, 0, MPI_TAG_MESH + 1, new NgMPI_Comm(comm));
        
			int maxidentnr = pp_data[0];
			var idents = GetIdentifications();
			for (int idnr = 1; idnr < maxidentnr + 1; idnr++)
			{
        
				idents.SetType(idnr, (Identifications.ID_TYPE)pp_data[idnr]);
			}
        
			int offset = 2 * maxidentnr + 1;
			for (int idnr = 1; idnr < maxidentnr + 1; idnr++)
			{
				int npairs = pp_data[maxidentnr + idnr];
				FlatArray<int> pairdata = new FlatArray<int>((uint)(2 * npairs), pp_data[offset]);
				offset += 2 * npairs;
			for (int k = 0; k < npairs; k++)
			{
			  PointIndex loc1 = glob2loc_vert_ht.Get(pairdata[2 * k]);
			  PointIndex loc2 = glob2loc_vert_ht.Get(pairdata[2 * k + 1]);
			  idents.Add(loc1, loc2, idnr);
			}
			}
        
			Array<int> dist_pnums = new Array<int>();
			MyMPI_Recv(ref dist_pnums, 0, MPI_TAG_MESH + 1, new NgMPI_Comm(comm));
        
			for (int hi = 0; hi < dist_pnums.Size(); hi += 3)
			{
			  paralleltop.SetDistantPNum(dist_pnums[hi + 1], dist_pnums[hi]); // , dist_pnums[hi+2]);
			}
        
			NgProfiler.StopTimer(timer_pts);
			*testout << "got " << numvert << " vertices" << "\n";
        
			{
			  Element el = new Element();
        
			  Array<int> elarray = new Array<int>();
			  MyMPI_Recv(ref elarray, 0, MPI_TAG_MESH + 2, new NgMPI_Comm(comm));
        
			  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer_els);
        
			  for (int ind = 0, elnum = 1; ind < elarray.Size(); elnum++)
			  {
			  paralleltop.SetLoc2Glob_VolEl(elnum, elarray[ind++]);
        
			  el.SetIndex(elarray[ind++]);
			  el.SetNP(elarray[ind++]);
        
			  for (int j = 0; j < el.GetNP(); j++)
			  {
				el[j] = glob2loc_vert_ht.Get(elarray[ind++]);
			  }
        
			  AddVolumeElement(el);
			  }
			}
        
			{
			  Array<double> fddata = new Array<double>();
			  MyMPI_Recv(ref fddata, 0, MPI_TAG_MESH + 3, new NgMPI_Comm(comm));
			  for (int i = 0; i < fddata.Size(); i += 6)
			  {
			  int faceind = AddFaceDescriptor(new FaceDescriptor((int)fddata[i], (int)fddata[i + 1], (int)fddata[i + 2], 0));
			  GetFaceDescriptor(faceind).SetBCProperty((int)fddata[i + 3]);
			  GetFaceDescriptor(faceind).domin_singular = fddata[i + 4];
			  GetFaceDescriptor(faceind).domout_singular = fddata[i + 5];
			  }
			}
        
			{
			  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer_sels);
			  Array<int> selbuf = new Array<int>();
        
			  MyMPI_Recv(ref selbuf, 0, MPI_TAG_MESH + 4, new NgMPI_Comm(comm));
        
			  int ii = 0;
			  int sel = 0;
        
			  int nlocsel = selbuf[ii++];
			  paralleltop.SetNSE(nlocsel);
        
			  while (ii < selbuf.Size() - 1)
			  {
			  int globsel = selbuf[ii++];
			  int faceind = selbuf[ii++];
			  //bool isghost = selbuf[ii++];
			  int nep = selbuf[ii++];
			  Element2d tri = new Element2d(nep);
			  tri.SetIndex(faceind);
			  for (int j = 1; j <= nep; j++)
			  {
				  tri.PNum(j) = glob2loc_vert_ht.Get(selbuf[ii++]);
				  tri.GeomInfoPi(j).trignum = selbuf[ii++];
			  }
			  paralleltop.SetLoc2Glob_SurfEl(sel + 1, globsel);
			  AddSurfaceElement(tri);
			  sel++;
			  }
			}
        
        
        
			{
			  Array<double> segmbuf = new Array<double>();
			  MyMPI_Recv(ref segmbuf, 0, MPI_TAG_MESH + 5, new NgMPI_Comm(comm));
        
			  Segment seg = new Segment();
			  int globsegi;
			  int ii = 0;
			  int segi = 1;
			  int nsegloc = segmbuf.Size() / 14;
			  paralleltop.SetNSegm(nsegloc);
        
			  while (ii < segmbuf.Size())
			  {
			  globsegi = (int)segmbuf[ii++];
			  seg.si = (int)segmbuf[ii++];
        
			  seg.pnums[0] = glob2loc_vert_ht.Get((int)segmbuf[ii++]);
			  seg.pnums[1] = glob2loc_vert_ht.Get((int)segmbuf[ii++]);
			  seg.geominfo[0].trignum = (int)segmbuf[ii++];
			  seg.geominfo[1].trignum = (int)segmbuf[ii++];
			  seg.surfnr1 = (int)segmbuf[ii++];
			  seg.surfnr2 = (int)segmbuf[ii++];
			  seg.edgenr = (int)segmbuf[ii++];
			  seg.epgeominfo[0].dist = segmbuf[ii++];
			  seg.epgeominfo[1].edgenr = (int)segmbuf[ii++];
			  seg.epgeominfo[1].dist = segmbuf[ii++];
        
			  seg.singedge_left = segmbuf[ii++];
			  seg.singedge_right = segmbuf[ii++];
        
			  seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;
        
			  seg.domin = seg.surfnr1;
			  seg.domout = seg.surfnr2;
			  if (seg.pnums[0] > 0 && seg.pnums[1] > 0)
			  {
				  paralleltop.SetLoc2Glob_Segm(segi, globsegi);
        
				  AddSegment(seg);
				  segi++;
			  }
        
			  }
			}
        
			/** Recv bc-names **/
			int[] nnames = {0, 0, 0, 0};
			MPI_Recv(nnames, 4, MPI_INT, 0, MPI_TAG_MESH + 6, comm, MPI_STATUS_IGNORE);
			materials.SetSize(nnames[0]);
			bcnames.SetSize(nnames[1]);
			cd2names.SetSize(nnames[2]);
			cd3names.SetSize(nnames[3]);
        
			int tot_nn = nnames[0] + nnames[1] + nnames[2] + nnames[3];
			Array<int> name_sizes = new Array<int>(tot_nn);
			MPI_Recv(name_sizes[0], tot_nn, MPI_INT, 0, MPI_TAG_MESH + 6, comm, MPI_STATUS_IGNORE);
			int tot_size = 0;
			for (int k = 0; k < tot_nn; k++)
			{
				tot_size += name_sizes[k];
			}
        
			Array<char> compiled_names = new Array<char>(tot_size);
			MPI_Recv((compiled_names[0]), tot_size, MPI_CHAR, 0, MPI_TAG_MESH + 6, comm, MPI_STATUS_IGNORE);
        
			tot_nn = tot_size = 0;
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto write_names = [&] (auto & array)
			var write_names = (array) =>
			{
			  for (int k = 0; k < array.Size(); k++)
			  {
			int s = name_sizes[tot_nn];
			array[k] = new string(compiled_names[tot_size], s);
			tot_nn++;
			tot_size += s;
			  }
			};
			write_names(materials);
			write_names(bcnames);
			write_names(cd2names);
			write_names(cd3names);
        
			MPI_Barrier(comm);
        
			int timerloc = NgProfiler.CreateTimer("Update local mesh");
			int timerloc2 = NgProfiler.CreateTimer("CalcSurfacesOfNode");
        
			NgProfiler.RegionTimer regloc = new NgProfiler.RegionTimer(timerloc);
			stringstream str = new stringstream();
			str << "p" << id << ": got " << GetNE() << " elements and " << GetNSE() << " surface elements";
			PrintMessage(2, str.str());
			// cout << str.str() << endl;
			// PrintMessage (2, "Got ", GetNE(), " elements and ", GetNSE(), " surface elements");
			// PrintMessage (2, "Got ", GetNSE(), " surface elements");
        
			NgProfiler.StartTimer(timerloc2);
        
			CalcSurfacesOfNode();
        
			NgProfiler.StopTimer(timerloc2);
        
			topology.Update();
			clusters.Update();
        
			// paralleltop -> UpdateCoarseGrid();
        
			SetNextMajorTimeStamp();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Distribute()
		  {
			NgMPI_Comm comm = GetCommunicator();
			int id = comm.Rank();
			int ntasks = comm.Size();
        
			if (id != 0 || ntasks == 1)
			{
				return;
			}
        
		#if METIS
			ParallelMetis();
		#else
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			{
			  this[ei].SetPartition(ntasks * ei / GetNE() + 1);
			}
		#endif
        
			/*
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			  *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
			for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
			  *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
			  */
        
			// MyMPI_SendCmd ("mesh");
			SendRecvMesh();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ParallelMetis()
		  {
			PrintMessage(3, "call metis 5 ...");
        
			int timer = NgProfiler.CreateTimer("Mesh::Partition");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer);
        
			idx_t ne = GetNE() + GetNSE() + GetNSeg();
			idx_t nn = GetNP();
        
			Array<idx_t> eptr = new Array<idx_t>();
			Array<idx_t> eind = new Array<idx_t>();
			for (int i = 0; i < GetNE(); i++)
			{
			eptr.Append(eind.Size());
			Element el = VolumeElement(i + 1);
			for (int j = 0; j < el.GetNP(); j++)
			{
			  eind.Append(el[j] - 1);
			}
			}
			for (int i = 0; i < GetNSE(); i++)
			{
			eptr.Append(eind.Size());
			Element2d el = SurfaceElement(i + 1);
			for (int j = 0; j < el.GetNP(); j++)
			{
			  eind.Append(el[j] - 1);
			}
			}
			for (int i = 0; i < GetNSeg(); i++)
			{
			eptr.Append(eind.Size());
			Segment el = LineSegment(i + 1);
			eind.Append(el[0] - 1);
			eind.Append(el[1] - 1);
			}
			eptr.Append(eind.Size());
			Array<idx_t> epart = new Array<idx_t>(ne);
			Array<idx_t> npart = new Array<idx_t>(nn);
        
			idxtype nparts = GetCommunicator().Size() - 1;
        
			vol_partition.SetSize(GetNE());
			surf_partition.SetSize(GetNSE());
			seg_partition.SetSize(GetNSeg());
			if (nparts == 1)
			{
				for (int i = 0; i < GetNE(); i++)
				{
				  // VolumeElement(i+1).SetPartition(1);
				  vol_partition[i] = 1;
				}
				for (int i = 0; i < GetNSE(); i++)
				{
				  // SurfaceElement(i+1).SetPartition(1);
				  surf_partition[i] = 1;
				}
				for (int i = 0; i < GetNSeg(); i++)
				{
				  // LineSegment(i+1).SetPartition(1);
				  seg_partition[i] = 1;
				}
			}
        
			else
        
			{
        
				idxtype edgecut = new idxtype();
        
				idxtype ncommon = 3;
				METIS_PartMeshDual(ne, nn, eptr[0], eind[0], null, null, ncommon, nparts, null, null, edgecut, epart[0], npart[0]);
        
				/*
				  METIS_PartMeshNodal (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &nparts,
				  NULL, NULL,
				  &edgecut, &epart[0], &npart[0]);
				*/
				PrintMessage(3, "metis complete");
				// cout << "done" << endl;
        
				for (int i = 0; i < GetNE(); i++)
				{
				  // VolumeElement(i+1).SetPartition(epart[i] + 1);
				  vol_partition[i] = epart[i] + 1;
				}
				for (int i = 0; i < GetNSE(); i++)
				{
				  // SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
				  surf_partition[i] = epart[i + GetNE()] + 1;
				}
				for (int i = 0; i < GetNSeg(); i++)
				{
				  // LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);
				  seg_partition[i] = epart[i + GetNE() + GetNSE()] + 1;
				}
			}
        
        
			// surface elements attached to volume elements
			Array<bool, PointIndex.BASE> boundarypoints = new Array<bool, PointIndex.BASE>(GetNP());
			boundarypoints = false;
        
			if (GetDimension() == 3)
			{
			  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			  {
			  Element2d el = this[sei];
			  for (int j = 0; j < el.GetNP(); j++)
			  {
				boundarypoints[el[j]] = true;
			  }
			  }
			}
			else
			{
			  for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
			  {
			  Segment seg = this[segi];
			  for (int j = 0; j < 2; j++)
			  {
				boundarypoints[seg[j]] = true;
			  }
			  }
			}
        
        
			// Build Pnt2Element table, boundary points only
			Array<int, PointIndex.BASE> cnt = new Array<int, PointIndex.BASE>(GetNP());
			cnt = 0;
        
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto loop_els_2d = [&](auto f)
			var loop_els_2d = (f UnnamedParameter) =>
			{
			  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			  {
			  Element2d el = this[sei];
			  for (int j = 0; j < el.GetNP(); j++)
			  {
				f(el[j], sei);
			  }
			  }
			};
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto loop_els_3d = [&](auto f)
			var loop_els_3d = (f UnnamedParameter) =>
			{
			  for (ElementIndex ei = 0; ei < GetNE(); ei++)
			  {
			  Element el = this[ei];
			  for (int j = 0; j < el.GetNP(); j++)
			  {
				f(el[j], ei);
			  }
			  }
			};
		//C++ TO C# CONVERTER TODO TASK: Lambda expressions cannot be assigned to 'var':
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: auto loop_els = [&](auto f)
			var loop_els = (f UnnamedParameter) =>
			{
			if (GetDimension() == 3)
			{
			  loop_els_3d(f);
			}
			else
			{
			  loop_els_2d(f);
			}
			};
        
        
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: loop_els([&](auto vertex, int index)
			loop_els((vertex UnnamedParameter, int index) =>
			{
			  if (boundarypoints[vertex])
			  {
				cnt[vertex]++;
			  }
			});
			TABLE<int, PointIndex.BASE> pnt2el = new TABLE<int, PointIndex.BASE>(cnt);
		//C++ TO C# CONVERTER NOTE: 'auto' variable declarations are not supported in C#:
		//ORIGINAL LINE: loop_els([&](auto vertex, int index)
			loop_els((vertex UnnamedParameter, int index) =>
			{
			  if (boundarypoints[vertex])
			  {
				pnt2el.Add(vertex, index);
			  }
			});
        
        
			if (GetDimension() == 3)
			{
			for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			{
				Element2d sel = this[sei];
				PointIndex pi1 = sel[0];
				// FlatArray<ElementIndex> els = pnt2el[pi1];
				FlatArray<int> els = pnt2el[pi1];
        
				// sel.SetPartition (-1);
					surf_partition[sei] = -1;
        
				for (int j = 0; j < els.Size(); j++)
				{
				Element el = this[new ElementIndex(els[j])];
        
				bool hasall = true;
        
				for (int k = 0; k < sel.GetNP(); k++)
				{
					bool haspi = false;
					for (int l = 0; l < el.GetNP(); l++)
					{
					  if (sel[k] == el[l])
					  {
					haspi = true;
					  }
					}
        
					if (!haspi)
					{
						hasall = false;
					}
				}
        
				if (hasall)
				{
					// sel.SetPartition (el.GetPartition());
							surf_partition[sei] = vol_partition[new ElementIndex(els[j])];
					break;
				}
				}
				// if (sel.GetPartition() == -1)
					if (surf_partition[sei] == -1)
					{
				  cerr << "no volume element found" << "\n";
					}
			}
        
        
			for (SegmentIndex si = 0; si < GetNSeg(); si++)
			{
				Segment sel = this[si];
				PointIndex pi1 = sel[0];
				FlatArray<int> els = pnt2el[pi1];
        
				// sel.SetPartition (-1);
					seg_partition[si] = -1;
        
				for (int j = 0; j < els.Size(); j++)
				{
				Element el = this[new ElementIndex(els[j])];
        
				bool[] haspi = {false, false, false, false, false, false, false, false, false}; // max surfnp
        
				for (int k = 0; k < 2; k++)
				{
				  for (int l = 0; l < el.GetNP(); l++)
				  {
					if (sel[k] == el[l])
					{
					  haspi[k] = true;
					}
				  }
				}
        
				bool hasall = true;
				for (int k = 0; k < sel.GetNP(); k++)
				{
				  if (!haspi[k])
				  {
					  hasall = false;
				  }
				}
        
				if (hasall)
				{
					// sel.SetPartition (el.GetPartition());
							seg_partition[si] = vol_partition[new ElementIndex(els[j])];
					break;
				}
				}
				// if (sel.GetPartition() == -1)
					if (seg_partition[si] == -1)
					{
				  cerr << "no volume element found" << "\n";
					}
			}
			}
			else
			{
			for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
			{
				Segment seg = this[segi];
				// seg.SetPartition(-1);
					seg_partition[segi] = -1;
				PointIndex pi1 = seg[0];
        
				FlatArray<int> sels = pnt2el[pi1];
				for (int j = 0; j < sels.Size(); j++)
				{
				SurfaceElementIndex sei = sels[j];
				Element2d se = this[sei];
				bool found = false;
				for (int l = 0; l < se.GetNP(); l++ && !found)
				{
				  found |= (se[l] == seg[1]);
				}
				if (found)
				{
				  // seg.SetPartition(se.GetPartition());
						  seg_partition[segi] = surf_partition[sei];
				  break;
				}
				}
        
				// if (seg.GetPartition() == -1) {
					if (seg_partition[segi] == -1)
					{
				  Console.Write("\n");
				  Console.Write("segi: ");
				  Console.Write(segi);
				  Console.Write("\n");
				  Console.Write("points: ");
				  Console.Write(seg[0]);
				  Console.Write(" ");
				  Console.Write(seg[1]);
				  Console.Write("\n");
				  Console.Write("surfels: ");
				  Console.Write("\n");
				  Console.Write(sels);
				  Console.Write("\n");
				  throw new Exception("no surface element found");
					}
			}
        
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Distribute(Array<int> volume_weights, Array<int> surface_weights, Array<int> segment_weights)
		  {
			NgMPI_Comm comm = GetCommunicator();
			int id = comm.Rank();
			int ntasks = comm.Size();
        
			if (id != 0 || ntasks == 1)
			{
				return;
			}
        
		#if METIS
			ParallelMetis(volume_weights, surface_weights, segment_weights);
		#else
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			{
			  this[ei].SetPartition(ntasks * ei / GetNE() + 1);
			}
		#endif
        
			/*
			for (ElementIndex ei = 0; ei < GetNE(); ei++)
			  *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
			for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
			  *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
			  */
        
			// MyMPI_SendCmd ("mesh");
			SendRecvMesh();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ParallelMetis(Array<int> volume_weights, Array<int> surface_weights, Array<int> segment_weights)
		  {
			PrintMessage(3, "call metis 5 with weights ...");
        
			// cout << "segment_weights " << segment_weights << endl;
			// cout << "surface_weights " << surface_weights << endl;
			// cout << "volume_weights " << volume_weights << endl;
        
			int timer = NgProfiler.CreateTimer("Mesh::Partition");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer);
        
			idx_t ne = GetNE() + GetNSE() + GetNSeg();
			idx_t nn = GetNP();
        
			Array<idx_t> eptr = new Array<idx_t>();
			Array<idx_t> eind = new Array<idx_t>();
			Array<idx_t> nwgt = new Array<idx_t>();
			for (int i = 0; i < GetNE(); i++)
			{
			eptr.Append(eind.Size());
        
			Element el = VolumeElement(i + 1);
        
			int ind = el.GetIndex();
			if (volume_weights.Size() < ind)
			{
				nwgt.Append(0);
			}
			else
			{
				nwgt.Append(volume_weights[ind - 1]);
			}
        
			for (int j = 0; j < el.GetNP(); j++)
			{
			  eind.Append(el[j] - 1);
			}
			}
			for (int i = 0; i < GetNSE(); i++)
			{
			eptr.Append(eind.Size());
			Element2d el = SurfaceElement(i + 1);
        
        
			int ind = el.GetIndex();
			ind = GetFaceDescriptor(ind).BCProperty();
			if (surface_weights.Size() < ind)
			{
				nwgt.Append(0);
			}
			else
			{
				nwgt.Append(surface_weights[ind - 1]);
			}
        
        
			for (int j = 0; j < el.GetNP(); j++)
			{
			  eind.Append(el[j] - 1);
			}
			}
			for (int i = 0; i < GetNSeg(); i++)
			{
			eptr.Append(eind.Size());
        
			Segment el = LineSegment(i + 1);
        
			int ind = el.si;
			if (segment_weights.Size() < ind)
			{
				nwgt.Append(0);
			}
			else
			{
				nwgt.Append(segment_weights[ind - 1]);
			}
        
			eind.Append(el[0]);
			eind.Append(el[1]);
			}
        
			eptr.Append(eind.Size());
			Array<idx_t> epart = new Array<idx_t>(ne);
			Array<idx_t> npart = new Array<idx_t>(nn);
        
			idxtype nparts = GetCommunicator().Size() - 1;
			vol_partition.SetSize(GetNE());
			surf_partition.SetSize(GetNSE());
			seg_partition.SetSize(GetNSeg());
        
			if (nparts == 1)
			{
				for (int i = 0; i < GetNE(); i++)
				{
				  // VolumeElement(i+1).SetPartition(1);
				  vol_partition[i] = 1;
				}
				for (int i = 0; i < GetNSE(); i++)
				{
				  // SurfaceElement(i+1).SetPartition(1);
				  surf_partition[i] = 1;
				}
				for (int i = 0; i < GetNSeg(); i++)
				{
				  // LineSegment(i+1).SetPartition(1);
				  seg_partition[i] = 1;
				}
				return;
			}
        
        
			idxtype edgecut = new idxtype();
        
        
			idxtype ncommon = 3;
			METIS_PartMeshDual(ne, nn, eptr[0], eind[0], nwgt[0], null, ncommon, nparts, null, null, edgecut, epart[0], npart[0]);
			/*
			METIS_PartMeshNodal (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &nparts,
					 NULL, NULL,
					 &edgecut, &epart[0], &npart[0]);
			*/
			PrintMessage(3, "metis complete");
			// cout << "done" << endl;
        
			for (int i = 0; i < GetNE(); i++)
			{
			  // VolumeElement(i+1).SetPartition(epart[i] + 1);
			  vol_partition[i] = epart[i] + 1;
			}
			for (int i = 0; i < GetNSE(); i++)
			{
			  // SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
			  surf_partition[i] = epart[i + GetNE()] + 1;
			}
			for (int i = 0; i < GetNSeg(); i++)
			{
			  // LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);
			  seg_partition[i] = epart[i + GetNE() + GetNSE()] + 1;
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ParallelMetis()
		  {
			int timer = NgProfiler.CreateTimer("Mesh::Partition");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer);
        
			PrintMessage(3, "Metis called");
        
			if (GetDimension() == 2)
			{
			PartDualHybridMesh2D(); // neloc );
			return;
			}
        
        
			idx_t ne = GetNE();
			idx_t nn = GetNP();
        
			if (ntasks <= 2 || ne <= 1)
			{
				if (ntasks == 1)
				{
					return;
				}
        
				for (int i = 1; i <= ne; i++)
				{
				  VolumeElement(i).SetPartition(1);
				}
        
				for (int i = 1; i <= GetNSE(); i++)
				{
				  SurfaceElement(i).SetPartition(1);
				}
        
				return;
			}
        
        
			bool uniform_els = true;
        
			ELEMENT_TYPE elementtype = ELEMENT_TYPE.TET;
			for (int el = 1; el <= GetNE(); el++)
			{
			  if (VolumeElement(el).GetType() != elementtype)
			  {
			  uniform_els = false;
			  break;
			  }
			}
        
        
			if (!uniform_els)
			{
			PartHybridMesh();
			}
			else
			{
        
			// uniform (TET) mesh,  JS
			int npe = VolumeElement(1).GetNP();
			Array<idxtype> elmnts = new Array<idxtype>(ne * npe);
        
			int etype;
			if (elementtype == ELEMENT_TYPE.TET)
			{
			  etype = 2;
			}
			else if (elementtype == ELEMENT_TYPE.HEX)
			{
			  etype = 3;
			}
        
        
			for (int i = 1; i <= ne; i++)
			{
			  for (int j = 1; j <= npe; j++)
			  {
				elmnts[(i - 1) * npe + (j - 1)] = VolumeElement(i).PNum(j) - 1;
			  }
			}
        
			int numflag = 0;
			int nparts = ntasks - 1;
			int ncommon = 3;
			int edgecut;
			Array<idxtype> epart = new Array<idxtype>(ne);
			Array<idxtype> npart = new Array<idxtype>(nn);
        
			//     if ( ntasks == 1 ) 
			//       {
			// 	(*this) = *mastermesh;
			// 	nparts = 4;	   
			// 	metis :: METIS_PartMeshDual (&ne, &nn, elmnts, &etype, &numflag, &nparts,
			// 				     &edgecut, epart, npart);
			// 	cout << "done" << endl;
        
			// 	cout << "edge-cut: " << edgecut << ", balance: " << metis :: ComputeElementBalance(ne, nparts, epart) << endl;
        
			// 	for (int i=1; i<=ne; i++)
			// 	  {
			// 	    mastermesh->VolumeElement(i).SetPartition(epart[i-1]);
			// 	  }
        
			// 	return;
			//       }
        
        
			int timermetis = NgProfiler.CreateTimer("Metis itself");
			NgProfiler.StartTimer(timermetis);
        
		//C++ TO C# CONVERTER TODO TASK: The cout 'flush' manipulator is not converted by C++ to C# Converter:
		//ORIGINAL LINE: cout << "call metis(4)_PartMeshDual ... " << flush;
		#if METIS4
			Console.Write("call metis(4)_PartMeshDual ... ");
			METIS_PartMeshDual(ne, nn, elmnts[0], etype, numflag, nparts, edgecut, epart[0], npart[0]);
		#else
			Console.Write("call metis(5)_PartMeshDual ... ");
			Console.Write("\n");
			// idx_t options[METIS_NOPTIONS];
        
			Array<idx_t> eptr = new Array<idx_t>(ne+1);
			for (int j = 0; j < ne+1; j++)
			{
			  eptr[j] = 4 * j;
			}
        
			METIS_PartMeshDual(ne, nn, eptr[0], elmnts[0], null, null, ncommon, nparts, null, null, edgecut, epart[0], npart[0]);
		#endif
        
			NgProfiler.StopTimer(timermetis);
        
			Console.Write("complete");
			Console.Write("\n");
		#if METIS4
			Console.Write("edge-cut: ");
			Console.Write(edgecut);
			Console.Write(", balance: ");
			Console.Write(ComputeElementBalance(ne, nparts, epart[0]));
			Console.Write("\n");
		#endif
        
			// partition numbering by metis : 0 ...  ntasks - 1
			// we want:                       1 ...  ntasks
			for (int i = 1; i <= ne; i++)
			{
			  VolumeElement(i).SetPartition(epart[i - 1] + 1);
			}
			}
        
        
			for (int sei = 1; sei <= GetNSE(); sei++)
			{
			int ei1;
			int ei2;
			GetTopology().GetSurface2VolumeElement(sei, ei1, ei2);
			Element2d sel = SurfaceElement(sei);
        
				for (int j = 0; j < 2; j++)
				{
					int ei = (j == 0) ? ei1 : ei2;
					if (ei > 0 && ei <= GetNE())
					{
				sel.SetPartition(VolumeElement(ei).GetPartition());
				break;
					}
				}
			}
        
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void PartHybridMesh()
		  {
		#if METIS
			int ne = GetNE();
        
			int nn = GetNP();
			int nedges = topology.GetNEdges();
        
			idxtype[] xadj;
			idxtype adjacency;
			idxtype v_weights = null;
			idxtype e_weights = null;
        
			int weightflag = 0;
			int numflag = 0;
			int nparts = ntasks - 1;
        
			int[] options = new int[5];
			options[0] = 0;
			int edgecut;
			idxtype[] part;
        
			xadj = Arrays.InitializeWithDefaultInstances<idxtype>(nn + 1);
			part = Arrays.InitializeWithDefaultInstances<idxtype>(nn);
        
			Array<int> cnt = new Array<int>(nn + 1);
			cnt = 0;
        
			for (int edge = 1; edge <= nedges; edge++)
			{
			int v1;
			int v2;
			topology.GetEdgeVertices(edge, v1, v2);
			cnt[v1 - 1]++;
			cnt[v2 - 1]++;
			}
        
			xadj[0] = 0;
			for (int n = 1; n <= nn; n++)
			{
			xadj[n] = idxtype(xadj[n - 1] + cnt[n - 1]);
			}
        
			adjacency = Arrays.InitializeWithDefaultInstances<idxtype>(xadj[nn]);
			cnt = 0;
        
			for (int edge = 1; edge <= nedges; edge++)
			{
			int v1;
			int v2;
			topology.GetEdgeVertices(edge, v1, v2);
			adjacency[xadj[v1 - 1] + cnt[v1 - 1]] = v2 - 1;
			adjacency[xadj[v2 - 1] + cnt[v2 - 1]] = v1 - 1;
			cnt[v1 - 1]++;
			cnt[v2 - 1]++;
			}
        
			for (int vert = 0; vert < nn; vert++)
			{
			FlatArray<idxtype> array = new FlatArray<idxtype>(cnt[vert], adjacency[xadj[vert]]);
			BubbleSort(array);
			}
        
		#if METIS4
			METIS_PartGraphKway(nn, xadj, adjacency, v_weights, e_weights, weightflag, numflag, nparts, options, edgecut, part);
		#else
			Console.Write("currently not supported (metis5), A");
			Console.Write("\n");
		#endif
        
			Array<int> nodesinpart = new Array<int>(ntasks);
			vol_partition.SetSize(ne);
			for (int el = 1; el <= ne; el++)
			{
			Element volel = VolumeElement(el);
			nodesinpart = 0;
        
        
			int el_np = volel.GetNP();
			int partition = 0;
			for (int i = 0; i < el_np; i++)
			{
			  nodesinpart[part[volel[i] - 1] + 1]++;
			}
        
			for (int i = 1; i < ntasks; i++)
			{
			  if (nodesinpart[i] > nodesinpart[partition])
			  {
				partition = i;
			  }
			}
        
			// volel.SetPartition(partition);
				vol_partition[el - 1] = partition;
			}
        
			xadj = null;
			part = null;
			adjacency = null;
		#else
			Console.Write("parthybridmesh not available");
			Console.Write("\n");
		#endif
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void PartDualHybridMesh()
		  {
		#if METIS
			int ne = GetNE();
        
			// int nn = GetNP();
			// int nedges = topology->GetNEdges();
			int nfaces = topology.GetNFaces();
        
			idxtype[] xadj;
			idxtype adjacency;
			idxtype v_weights = null;
			idxtype e_weights = null;
        
			int weightflag = 0;
			// int numflag = 0;
			int nparts = ntasks - 1;
        
			int[] options = new int[5];
			options[0] = 0;
			int edgecut;
			idxtype[] part;
        
			Array<int, 0> facevolels1 = new Array<int, 0>(nfaces);
			Array<int, 0> facevolels2 = new Array<int, 0>(nfaces);
			facevolels1 = -1;
			facevolels2 = -1;
        
			Array<int, 0> elfaces = new Array<int, 0>();
			xadj = Arrays.InitializeWithDefaultInstances<idxtype>(ne+1);
			part = Arrays.InitializeWithDefaultInstances<idxtype>(ne);
        
			Array<int, 0> cnt = new Array<int, 0>(ne+1);
			cnt = 0;
        
			for (int el = 1; el <= ne; el++)
			{
			Element volel = VolumeElement(el);
			topology.GetElementFaces(el, elfaces);
			for (int i = 0; i < elfaces.Size(); i++)
			{
				if (facevolels1[elfaces[i] - 1] == -1)
				{
				  facevolels1[elfaces[i] - 1] = el;
				}
				else
				{
				facevolels2[elfaces[i] - 1] = el;
				cnt[facevolels1[elfaces[i] - 1] - 1]++;
				cnt[facevolels2[elfaces[i] - 1] - 1]++;
				}
			}
			}
        
			xadj[0] = 0;
			for (int n = 1; n <= ne; n++)
			{
			xadj[n] = idxtype(xadj[n - 1] + cnt[n - 1]);
			}
        
			adjacency = Arrays.InitializeWithDefaultInstances<idxtype>(xadj[ne]);
			cnt = 0;
        
			for (int face = 1; face <= nfaces; face++)
			{
			int e1;
			int e2;
			e1 = facevolels1[face-1];
			e2 = facevolels2[face-1];
			if (e2 == -1)
			{
				continue;
			}
			adjacency[xadj[e1 - 1] + cnt[e1 - 1]] = e2 - 1;
			adjacency[xadj[e2 - 1] + cnt[e2 - 1]] = e1 - 1;
			cnt[e1 - 1]++;
			cnt[e2 - 1]++;
			}
        
			for (int el = 0; el < ne; el++)
			{
			FlatArray<idxtype> array = new FlatArray<idxtype>(cnt[el], adjacency[xadj[el]]);
			BubbleSort(array);
			}
        
			int timermetis = NgProfiler.CreateTimer("Metis itself");
			NgProfiler.StartTimer(timermetis);
        
		#if METIS4
			METIS_PartGraphKway(ne, xadj, adjacency, v_weights, e_weights, weightflag, numflag, nparts, options, edgecut, part);
		#else
			Console.Write("currently not supported (metis5), B");
			Console.Write("\n");
		#endif
        
        
			NgProfiler.StopTimer(timermetis);
        
			Array<int> nodesinpart = new Array<int>(ntasks);
        
			vol_partition.SetSize(ne);
			for (int el = 1; el <= ne; el++)
			{
			// Element & volel = VolumeElement(el);
			nodesinpart = 0;
        
			// VolumeElement(el).SetPartition(part[el-1 ] + 1);
			vol_partition[el - 1] = part[el - 1] + 1;
			}
        
			/*    
			for ( int i=1; i<=ne; i++)
			  {
			neloc[ VolumeElement(i).GetPartition() ] ++;
			  }
			*/
        
			xadj = null;
			part = null;
			adjacency = null;
		#else
			Console.Write("partdualmesh not available");
			Console.Write("\n");
		#endif
        
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void PartDualHybridMesh2D()
		  {
		#if METIS
			idxtype ne = GetNSE();
			int nv = GetNV();
        
			Array<idxtype> xadj = new Array<idxtype>(ne+1);
			Array<idxtype> adjacency = new Array<idxtype>(ne * 4);
        
			// first, build the vertex 2 element table:
			Array<int, PointIndex.BASE> cnt = new Array<int, PointIndex.BASE>(nv);
			cnt = 0;
			for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			{
			  for (int j = 0; j < this[sei].GetNP(); j++)
			  {
			cnt[this[sei][j]]++;
			  }
			}
        
			TABLE<SurfaceElementIndex, PointIndex.BASE> vert2els = new TABLE<SurfaceElementIndex, PointIndex.BASE>(cnt);
			for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			{
			  for (int j = 0; j < this[sei].GetNP(); j++)
			  {
			vert2els.Add(this[sei][j], sei);
			  }
			}
        
        
			// find all neighbour elements
			int cntnb = 0;
			Array<int> marks = new Array<int>(ne); // to visit each neighbour just once
			marks = -1;
			for (SurfaceElementIndex sei = 0; sei < ne; sei++)
			{
			xadj[sei] = cntnb;
			for (int j = 0; j < this[sei].GetNP(); j++)
			{
				PointIndex vnr = this[sei][j];
        
				// all elements with at least one common vertex
				for (int k = 0; k < vert2els[vnr].Size(); k++)
				{
				SurfaceElementIndex sei2 = vert2els[vnr][k];
				if (sei == sei2)
				{
					continue;
				}
				if (marks[sei2] == sei)
				{
					continue;
				}
        
				// neighbour, if two common vertices
				int common = 0;
				for (int m1 = 0; m1 < this[sei].GetNP(); m1++)
				{
				  for (int m2 = 0; m2 < this[sei2].GetNP(); m2++)
				  {
					if (this[sei][m1] == this[sei2][m2])
					{
					  common++;
					}
				  }
				}
        
				if (common >= 2)
				{
					marks[sei2] = sei; // mark as visited
					adjacency[cntnb++] = sei2;
				}
				}
			}
			}
			xadj[ne] = cntnb;
        
			idxtype v_weights = null;
			idxtype e_weights = null;
        
			idxtype weightflag = 0;
			// int numflag = 0;
			idxtype nparts = ntasks - 1;
        
			idxtype edgecut = new idxtype();
			Array<idxtype> part = new Array<idxtype>(ne);
        
			for (int el = 0; el < ne; el++)
			{
			  BubbleSort(adjacency.Range(xadj[el], xadj[el + 1]));
			}
        
		#if METIS4
			int[] options = new int[5];
			options[0] = 0;
			METIS_PartGraphKway(ne, xadj[0], adjacency[0], v_weights, e_weights, weightflag, numflag, nparts, options, edgecut, part[0]);
		#else
			idx_t ncon = 1;
			METIS_PartGraphKway(ne, ncon, xadj[0], adjacency[0], v_weights, null, e_weights, nparts, null, null, null, edgecut, part[0]);
		#endif
        
        
			surf_partition.SetSize(ne);
			for (SurfaceElementIndex sei = 0; sei < ne; sei++)
			{
			  // (*this) [sei].SetPartition (part[sei]+1);
			  surf_partition[sei] = part[sei] + 1;
			}
		#else
			Console.Write("partdualmesh not available");
			Console.Write("\n");
		#endif
        
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void ImproveMesh(CSG eometry geometry, OPTIMIZEGOAL goal)
		{
		  int i;
		  int eli;
		  int j;
		  int typ = 1;
        
		  if (geometry == null || geometry.GetNSurf() == 0)
		  {
			  ImproveMesh(goal);
			  return;
		  }
        
		  string savetask = multithread.task;
		  multithread.task = "Smooth Mesh";
        
        
		  TABLE<int> surfelementsonpoint = new TABLE<int>(points.Size());
		  Vector x = new Vector(3);
		  Vector xsurf = new Vector(2);
		  Vector xedge = new Vector(1);
		  int surf;
		  int surf1;
		  int surf2;
		  int surf3;
        
		  int uselocalh = mparam.uselocalh;
        
		  (*testout) << setprecision(8);
		  (*testout) << "Improve Mesh" << "\n";
		  PrintMessage(3, "ImproveMesh");
		  //  (*mycout) << "Vol = " << CalcVolume (points, volelements) << endl;
        
        
		  for (i = 1; i <= surfelements.Size(); i++)
		  {
			for (j = 1; j <= 3; j++)
			{
			  surfelementsonpoint.Add1(surfelements.Get(i).PNum(j), i);
			}
		  }
        
        
		  PointFunction pf;
		  if (typ == 1)
		  {
			pf = new PointFunction(points, volelements);
		  }
		  else
		  {
			pf = new CheapPointFunction(points, volelements);
		  }
        
		  //  pf->SetLocalH (h);
        
		  Opti3FreeMinFunction freeminf = new Opti3FreeMinFunction(pf);
		  Opti3SurfaceMinFunction surfminf = new Opti3SurfaceMinFunction(pf);
		  Opti3EdgeMinFunction edgeminf = new Opti3EdgeMinFunction(pf);
        
		  OptiParameters par = new OptiParameters();
		  par.maxit_linsearch = 20;
		  par.maxit_bfgs = 20;
        
		  int printmod = 1;
		  char printdot = '.';
		  if (points.Size() > 1000)
		  {
			  printmod = 10;
			  printdot = '+';
		  }
		  if (points.Size() > 10000)
		  {
			  printmod = 100;
			  printdot = '*';
		  }
        
		  for (i = 1; i <= points.Size(); i++)
		  {
			  //      if (ptyps.Get(i) == FIXEDPOINT) continue;
			  if (ptyps.Get(i) != POINTTYPE.INNERPOINT)
			  {
				  continue;
			  }
        
			  if (multithread.terminate)
			  {
			throw new Exception("Meshing stopped");
			  }
			  /*
			  if (multithread.terminate)
			break;
			  */
			  multithread.percent = 100.0 * i / points.Size();
        
			  /*
			  if (points.Size() < 1000)
			PrintDot ();
			  else
			if (i % 10 == 0)
			  PrintDot ('+');
			  */
			  if (i % printmod == 0)
			  {
				  PrintDot(printdot);
			  }
        
			  //    (*testout) << "Now point " << i << "\n";
			  //    (*testout) << "Old: " << points.Get(i) << "\n";
        
			  pf.SetPointIndex(i);
        
			  {
			  //      if (uselocalh)
			double lh = GetH(points.Get(i));
			pf.SetLocalH(GetH(points.Get(i)));
			par.typx = lh / 10;
			//	  (*testout) << "lh(" << points.Get(i) << ") = " << lh << "\n";
			  }
        
			  surf1 = surf2 = surf3 = 0;
        
			  for (j = 1; j <= surfelementsonpoint.EntrySize(i); j++)
			  {
			  eli = surfelementsonpoint.Get(i, j);
			  int surfi = surfelements.Get(eli).GetIndex();
        
			  if (surfi != 0)
			  {
				  surf = GetFaceDescriptor(surfi).SurfNr();
        
				  if (surf1 == 0)
				  {
				surf1 = surf;
				  }
				  else if (surf1 != surf)
				  {
				  if (surf2 == 0)
				  {
					surf2 = surf;
				  }
				  else if (surf2 != surf)
				  {
					surf3 = surf;
				  }
				  }
			  }
			  else
			  {
				  surf1 = surf2 = surf3 = 1; // simulates corner point
			  }
			  }
        
        
			  if (surf2 != 0 && surf3 == 0)
			  {
			  //      (*testout) << "On Edge" << "\n";
			  /*
			    xedge = 0;
			    edgeminf.SetPoint (geometry.GetSurface(surf1),
			    geometry.GetSurface(surf2), 
			    points.Elem(i));
			    BFGS (xedge, edgeminf, par);
		
			    edgeminf.CalcNewPoint (xedge, points.Elem(i));
			  */
			  }
        
			  if (surf1 != 0 && surf2 == 0)
			  {
			  //      (*testout) << "In Surface" << "\n";
			  /*
			    xsurf = 0;
			    surfminf.SetPoint (geometry.GetSurface(surf1),
			    points.Get(i));
			    BFGS (xsurf, surfminf, par);
		
			    surfminf.CalcNewPoint (xsurf, points.Elem(i));
			  */
			  }
        
			  if (surf1 == 0)
			  {
			  //      (*testout) << "In Volume" << "\n";
			  x = 0;
			  freeminf.SetPoint(points.Elem(i));
			  //	  par.typx = 
			  BFGS(x, freeminf, par);
        
			  points.Elem(i).X() += x.Get(1);
			  points.Elem(i).Y() += x.Get(2);
			  points.Elem(i).Z() += x.Get(3);
			  }
        
			  //    (*testout) << "New Point: " << points.Elem(i) << "\n" << "\n";
        
		  }
		  PrintDot('\n');
		  //  (*mycout) << "Vol = " << CalcVolume (points, volelements) << endl;
        
		  multithread.task = savetask;
        
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void ImproveMesh(MeshingParameters mp, OPTIMIZEGOAL goal)
		{
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//  static Timer t("Mesh::ImproveMesh");
		  RegionTimer reg = new RegionTimer(ImproveMesh_t);
        
		  int typ = 1;
        
		  (*testout) << "Improve Mesh" << "\n";
		  PrintMessage(3, "ImproveMesh");
        
		  int np = GetNP();
		  int ne = GetNE();
        
        
		  Array<double,PointIndex.BASE> perrs = new Array<double,PointIndex.BASE>(np);
		  perrs = 1.0;
        
		  double bad1 = 0;
		  double badmax = 0;
        
		  if (goal == OPTIMIZEGOAL.OPT_QUALITY)
		  {
			  for (int i = 1; i <= ne; i++)
			  {
			  Element el = VolumeElement(i);
			  if (el.GetType() != ELEMENT_TYPE.TET)
			  {
				continue;
			  }
        
			  double hbad = CalcBad(points, el, 0, mp);
			  for (int j = 0; j < 4; j++)
			  {
				perrs[el[j]] += hbad;
			  }
        
			  bad1 += hbad;
			  }
        
			  for (int i = perrs.Begin(); i < perrs.End(); i++)
			  {
			if (perrs[i] > badmax)
			{
			  badmax = perrs[i];
			}
			  }
			  badmax = 0;
		  }
        
		  if (goal == OPTIMIZEGOAL.OPT_QUALITY)
		  {
			  bad1 = CalcTotalBad(points, volelements, mp);
			  (*testout) << "Total badness = " << bad1 << "\n";
			  PrintMessage(5, "Total badness = ", bad1);
		  }
        
		  Vector x = new Vector(3);
        
		  (*testout) << setprecision(8);
        
		  //int uselocalh = mparam.uselocalh;
        
        
		  PointFunction pf;
        
		  if (typ == 1)
		  {
			pf = new PointFunction(points, volelements, mp);
		  }
		  else
		  {
			pf = new CheapPointFunction(points, volelements, mp);
		  }
        
		  //  pf->SetLocalH (h);
        
		  Opti3FreeMinFunction freeminf = new Opti3FreeMinFunction(pf);
        
		  OptiParameters par = new OptiParameters();
		  par.maxit_linsearch = 20;
		  par.maxit_bfgs = 20;
        
		  Array<double, PointIndex.BASE> pointh = new Array<double, PointIndex.BASE>(points.Size());
        
		  if (lochfunc)
		  {
			  foreach (PointIndex pi in points.Range())
			  {
			pointh[pi] = GetH(points[pi]);
			  }
		  }
		  else
		  {
			  pointh = 0;
			  foreach (Element el in VolumeElements())
			  {
			  double h = ngsimd.GlobalMembers.pow(el.Volume(points), 1.0 / 3.0);
				  foreach (PointIndex pi in el.PNums())
				  {
				if (h > pointh[pi])
				{
					  pointh[pi] = h;
				}
				  }
			  }
		  }
        
        
		  int printmod = 1;
		  char printdot = '.';
		  if (points.Size() > 1000)
		  {
			  printmod = 10;
			  printdot = '+';
		  }
		  if (points.Size() > 10000)
		  {
			  printmod = 100;
			  printdot = '*';
		  }
        
        
		  string savetask = multithread.task;
		  multithread.task = "Smooth Mesh";
        
		  foreach (PointIndex pi in points.Range())
		  {
			if (this[pi].Type() == POINTTYPE.INNERPOINT && perrs[pi] > 0.01 * badmax)
			{
			if (multithread.terminate)
			{
			  throw new Exception("Meshing stopped");
			}
        
			multithread.percent = 100.0 * (pi + 1 - PointIndex.BASE) / points.Size();
        
				if ((pi + 1 - PointIndex.BASE) % printmod == 0 != null)
				{
					PrintDot(printdot);
				}
        
			double lh = pointh[pi];
			pf.SetLocalH(lh);
			par.typx = lh;
        
			freeminf.SetPoint(points[pi]);
			pf.SetPointIndex(new netgen.PointIndex(pi));
        
			x = 0;
			int pok;
			pok = freeminf.Func(x) < 1e10;
        
			if (pok == 0)
			{
				pok = pf.MovePointToInner();
        
				freeminf.SetPoint(points[pi]);
				pf.SetPointIndex(new netgen.PointIndex(pi));
			}
        
			if (pok != 0)
			{
					//*testout << "start BFGS, pok" << endl;
				BFGS(x, freeminf, par);
					//*testout << "BFGS complete, pok" << endl;
				points[pi](0) += x(0);
				points[pi](1) += x(1);
				points[pi](2) += x(2);
			}
			}
		  }
		  PrintDot('\n');
        
		  if (pf != null)
		  {
			  pf.Dispose();
		  }
        
		  multithread.task = savetask;
        
		  if (goal == OPTIMIZEGOAL.OPT_QUALITY)
		  {
			  bad1 = CalcTotalBad(points, volelements, mp);
			  (*testout) << "Total badness = " << bad1 << "\n";
			  PrintMessage(5, "Total badness = ", bad1);
		  }
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void ImproveMeshJacobian(MeshingParameters mp, OPTIMIZEGOAL goal, BitArray usepoint)
		{
		  int i;
		  int j;
        
		  (*testout) << "Improve Mesh Jacobian" << "\n";
		  PrintMessage(3, "ImproveMesh Jacobian");
        
		  int np = GetNP();
		  int ne = GetNE();
        
        
		  Vector x = new Vector(3);
        
		  (*testout) << setprecision(8);
        
		  JacobianPointFunction pf = new JacobianPointFunction(points, volelements);
        
        
		  OptiParameters par = new OptiParameters();
		  par.maxit_linsearch = 20;
		  par.maxit_bfgs = 20;
        
		  BitArray badnodes = new BitArray(np);
		  badnodes.Clear();
        
		  for (i = 1; i <= ne; i++)
		  {
			  Element el = VolumeElement(i);
			  double bad = el.CalcJacobianBadness(Points());
			  if (bad > 1)
			  {
			for (j = 1; j <= el.GetNP(); j++)
			{
			  badnodes.Set(el.PNum(j));
			}
			  }
		  }
        
		  Array<double, PointIndex.BASE> pointh = new Array<double, PointIndex.BASE>(points.Size());
        
		  if (lochfunc)
		  {
			  for (i = 1; i <= points.Size(); i++)
			  {
			pointh[i] = GetH(points.Get(i));
			  }
		  }
		  else
		  {
			  pointh = 0;
			  for (i = 0; i < GetNE(); i++)
			  {
			  Element el = VolumeElement(i + 1);
			  double h = ngsimd.GlobalMembers.pow(el.Volume(points), 1.0 / 3.0);
			  for (j = 1; j <= el.GetNV(); j++)
			  {
				if (h > pointh[el.PNum(j)])
				{
				  pointh[el.PNum(j)] = h;
				}
			  }
			  }
		  }
        
        
        
		  string savetask = multithread.task;
		  multithread.task = "Smooth Mesh Jacobian";
        
		  for (PointIndex pi = points.Begin(); i < points.End(); pi++)
		  {
			  if (this[pi].Type() != POINTTYPE.INNERPOINT)
			  {
			continue;
			  }
        
			  if (usepoint != null && !usepoint.Test(i))
			  {
			continue;
			  }
        
			  //(*testout) << "improvejac, p = " << i << endl;
        
			  if (goal == OPTIMIZEGOAL.OPT_WORSTCASE && !badnodes.Test(i))
			  {
			continue;
			  }
			  //	(*testout) << "smooth p " << i << endl;
        
			  /*
			if (multithread.terminate)
			break;
			  */
			  if (multithread.terminate)
			  {
			throw new Exception("Meshing stopped");
			  }
        
			  multithread.percent = 100.0 * i / points.Size();
        
			  if (points.Size() < 1000)
			  {
			PrintDot();
			  }
			  else
			  {
			if (i % 10 == 0)
			{
			  PrintDot('+');
			}
			  }
        
			  double lh = pointh[i];
			  par.typx = lh;
        
			  pf.SetPointIndex(new netgen.PointIndex(pi));
        
			  x = 0;
			  int pok = (pf.Func(x) < 1e10);
        
			  if (pok != 0)
			  {
				  //*testout << "start BFGS, Jacobian" << endl;
			  BFGS(x, pf, par);
				  //*testout << "end BFGS, Jacobian" << endl;
			  points.Elem(i)(0) += x(0);
			  points.Elem(i)(1) += x(1);
			  points.Elem(i)(2) += x(2);
			  }
			  else
			  {
			  Console.Write("el not ok");
			  Console.Write("\n");
			  }
		  }
		  PrintDot('\n');
        
        
		  multithread.task = savetask;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void ImproveMeshJacobianOnSurface(MeshingParameters mp, BitArray usepoint, Array< Vec < 3> > nv, OPTIMIZEGOAL goal, Array< Array<int,PointIndex.BASE> > idmaps)
		{
		  int i;
		  int j;
        
		  (*testout) << "Improve Mesh Jacobian" << "\n";
		  PrintMessage(3, "ImproveMesh Jacobian");
        
		  int np = GetNP();
		  int ne = GetNE();
        
        
		  Vector x = new Vector(3);
        
		  testout.precision(8);
        
		  JacobianPointFunction pf = new JacobianPointFunction(points, volelements);
        
		  Array< Array<int,PointIndex.BASE> > locidmaps = new Array< Array<int,PointIndex.BASE> >();
		  Array< Array<int,PointIndex.BASE> > used_idmaps;
        
		  if (idmaps != null)
		  {
			used_idmaps = idmaps;
		  }
		  else
		  {
			  used_idmaps = locidmaps;
        
			  for (i = 1; i <= GetIdentifications().GetMaxNr(); i++)
			  {
			  if (GetIdentifications().GetType(i) == Identifications.PERIODIC)
			  {
				  locidmaps.Append(new Array<int,PointIndex.BASE>());
				  GetIdentifications().GetMap(i,*locidmaps.Last(),true);
			  }
			  }
		  }
        
        
		  bool usesum = (used_idmaps.Size() > 0);
		  MinFunctionSum pf_sum = new MinFunctionSum();
        
		  JacobianPointFunction pf2ptr = null;
		  if (usesum)
		  {
			  pf2ptr = new JacobianPointFunction(points, volelements);
			  pf_sum.AddFunction(pf);
			  pf_sum.AddFunction(pf2ptr);
		  }
        
        
		  OptiParameters par = new OptiParameters();
		  par.maxit_linsearch = 20;
		  par.maxit_bfgs = 20;
        
		  BitArray badnodes = new BitArray(np);
		  badnodes.Clear();
        
		  for (i = 1; i <= ne; i++)
		  {
			  Element el = VolumeElement(i);
			  double bad = el.CalcJacobianBadness(Points());
			  if (bad > 1)
			  {
			for (j = 1; j <= el.GetNP(); j++)
			{
			  badnodes.Set(el.PNum(j));
			}
			  }
		  }
        
		  Array<double, PointIndex.BASE> pointh = new Array<double, PointIndex.BASE>(points.Size());
        
		  if (lochfunc)
		  {
			  for (i = 1; i <= points.Size(); i++)
			  {
			pointh[i] = GetH(points.Get(i));
			  }
		  }
		  else
		  {
			  pointh = 0;
			  for (i = 0; i < GetNE(); i++)
			  {
			  Element el = VolumeElement(i + 1);
			  double h = ngsimd.GlobalMembers.pow(el.Volume(points), 1.0 / 3.0);
			  for (j = 1; j <= el.GetNV(); j++)
			  {
				if (h > pointh[el.PNum(j)])
				{
				  pointh[el.PNum(j)] = h;
				}
			  }
			  }
		  }
        
        
		  string savetask = multithread.task;
		  multithread.task = "Smooth Mesh Jacobian";
        
		  for (PointIndex pi = points.Begin(); pi <= points.End(); pi++)
		  {
			if (usepoint.Test(i))
			{
			//(*testout) << "improvejac, p = " << i << endl;
        
			if (goal == OPTIMIZEGOAL.OPT_WORSTCASE && !badnodes.Test(i))
			{
			  continue;
			}
			//	(*testout) << "smooth p " << i << endl;
        
			/*
			if (multithread.terminate)
			  break;
			*/
			if (multithread.terminate)
			{
			  throw new Exception("Meshing stopped");
			}
        
			multithread.percent = 100.0 * i / points.Size();
        
			if (points.Size() < 1000)
			{
			  PrintDot();
			}
			else
			{
			  if (i % 10 == 0)
			  {
				PrintDot('+');
			  }
			}
        
			double lh = pointh[i]; //GetH(points.Get(i));
			par.typx = lh;
        
			pf.SetPointIndex(new netgen.PointIndex(pi));
        
			PointIndex brother = new PointIndex(-1);
			if (usesum)
			{
				for (j = 0; brother == -1 && j < used_idmaps.Size(); j++)
				{
				if (i < used_idmaps[j].Size() + PointIndex.BASE)
				{
					brother.CopyFrom((* used_idmaps[j])[i]);
					if (brother == i || brother == 0)
					{
					  brother = -1;
					}
				}
				}
				if (brother >= i)
				{
				pf2ptr.SetPointIndex(new netgen.PointIndex(brother));
				pf2ptr.SetNV(nv[brother - 1]);
				}
			}
        
			if (usesum && brother < i)
			{
			  continue;
			}
        
			//pf.UnSetNV(); x = 0;
			//(*testout) << "before " << pf.Func(x);
        
			pf.SetNV(nv[i - 1]);
        
			x = 0;
			int pok = (brother == -1) ? (pf.Func(x) < 1e10) : (pf_sum.Func(x) < 1e10);
        
			if (pok != 0)
			{
        
				if (brother == -1)
				{
				  BFGS(x, pf, par);
				}
				else
				{
				  BFGS(x, pf_sum, par);
				}
        
        
				for (j = 0; j < 3; j++)
				{
				  points.Elem(i)(j) += x(j); // - scal*nv[i-1].X(j);
				}
        
				if (brother != -1)
				{
				  for (j = 0; j < 3; j++)
				  {
				points.Elem(brother)(j) += x(j); // - scal*nv[brother-1].X(j);
				  }
				}
        
        
			}
			else
			{
				Console.Write("el not ok");
				Console.Write("\n");
				(*testout) << "el not ok" << "\n" << "   func " << ((brother == -1) ? pf.Func(x) : pf_sum.Func(x)) << "\n";
				if (brother != -1)
				{
				  (*testout) << "   func1 " << pf.Func(x) << "\n" << "   func2 " << pf2ptr.Func(x) << "\n";
				}
			}
			}
		  }
        
		  PrintDot('\n');
        
		  if (pf2ptr != null)
		  {
			  pf2ptr.Dispose();
		  }
		  for (i = 0; i < locidmaps.Size(); i++)
		  {
			locidmaps[i] = null;
		  }
        
		  multithread.task = savetask;
		}
	}
}