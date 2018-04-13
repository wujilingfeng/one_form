#pragma once
#pragma once
#include<iostream>
#include <list>
#include<nanogui\screen.h>
//#include<nanogui\formhelper.h>
#include<igl\readOFF.h>
#include<igl\viewer\Viewer.h>
#include<igl\jet.h>
//#include<Eigen\Core>
//#include <Eigen\Eigenvalues>
//#include<OpenMesh\Core\IO\IOManager.hh>
#include<OpenMesh\Core\IO\MeshIO.hh>
#include<OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
#include<OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
#include <SymEigsSolver.h> 
#include<GenEigsSolver.h>
#include<MatOp/SparseGenMatProd.h>
using namespace Spectra;
using namespace std;
using namespace OpenMesh;
#ifndef N
#define N 13
#endif
struct MyTraits:OpenMesh::DefaultTraits
{
	VertexTraits{
		int layer = -1;
	float weight = -1;
	int key_ = 0;
	float f_ = 0;
	float area = 0;
	};
	EdgeTraits{
		float length = -1;
	bool sharp_ = false;
	};

};
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PMyMesh;

OpenMesh::HalfedgeHandle halfedgeh(VertexHandle v1, VertexHandle v2,PMyMesh &mesh)
{
	HalfedgeHandle he = mesh.halfedge_handle(v1), he1 = he;
	//mesh.edge
	if (mesh.from_vertex_handle(he1) != v1)
	{
		cout << "错误：：" << endl;
	}
	do
	{
		he1 = mesh.ccw_rotated_halfedge_handle(he1);
		//he1 = (he1);
	} while (he1 != he && v2 != mesh.to_vertex_handle(he1));
	return he1;
}
//typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PMyMesh;

class Mos
{
public:
	Mos(PMyMesh &mesh_, igl::viewer::Viewer &viewer_):mesh(mesh_),viewer(viewer_)
	{
		printf_s("create mos\r\n");
	}
	~Mos()
	{}
	void init_A()
	{
		int num_v = mesh.n_vertices();
		A.resize(num_v, num_v);
		A.setZero();
		for (PMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			//temp = 0;
			for (PMyMesh::VertexVertexIter vvter = mesh.vv_iter(*vter); vvter.is_valid(); vvter++)
			{
				//temp++;
				A.coeffRef((*vter).idx(), (*vvter).idx()) = 1;
				A.coeffRef((*vvter).idx(), (*vter).idx()) = 1;

			}
			int tempv = (mesh.valence((*vter)));
			//cout << "temp  " << temp << "  valence  " << mesh.valence((*vter)) << endl;
			A.coeffRef((*vter).idx(), (*vter).idx()) = -tempv;
			//cout << "temp " << temp <<"id  "<<(*vter).idx()<< endl;
		}
		
	}
	void set_f()
	{
		init_A();
		DenseSymMatProd<double> op(A);
		SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 2, 100);
		eigs.init();
		int nconv = eigs.compute();
		Eigen::VectorXd evalues;
		if (eigs.info() == SUCCESSFUL)
		{
			evalues = eigs.eigenvalues();}
		cout << evalues << endl;
		Eigen::MatrixXd eigenvectors;
		eigenvectors = eigs.eigenvectors(3);
		//cout << eigenvectors.col(1) << endl;
		//set_f
		for (PMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			mesh.data(vter).f_ = (1000 * eigenvectors.coeff((*vter).idx(), 1));
			//mesh.data(*vter).set_f((*vter).idx());
		}
	}
	void init_sharp()
	{

		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			mesh.data(*eter).sharp_ = (false);
		}
	}
	void init_edge_length()
	{
		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			mesh.data(eter).length = mesh.calc_edge_length(*eter);
			printf_s("%f \r\n", mesh.data(eter).length);
		}
	}
	void find_nextlayer(VertexHandle v)
	{
		VertexHandle v1;
		int layer = mesh.data(v).layer;
		float weight = mesh.data(v).weight;
		for (PMyMesh::VertexIHalfedgeIter vheter = mesh.vih_begin(v); vheter.is_valid(); vheter++)
		{
			v1 = mesh.from_vertex_handle(*vheter);
			if (mesh.data(v1).layer == -1)
			{
				mesh.data(v1).layer = layer + 1;
				mesh.data(v1).weight = mesh.data(mesh.edge_handle(*vheter)).length + weight;
				mylst1.push_back(v1);

			}
			else if (mesh.data(v1).layer == (layer + 1))
			{
				float temp = mesh.data(mesh.edge_handle(*vheter)).length + weight;
				if (temp<mesh.data(v1).weight)
				{
					mesh.data(v1).weight = temp;
				}
			}
			else if (mesh.data(v1).layer > layer+1)
			{
				printf_s("错误\r\n");
				return;
			}
		}

	}
	void set_destination(VertexHandle v)
	{
		mesh.data(v).layer = 0;
		mesh.data(v).weight = 0;
		mylst0.push_back(v);
		while (!mylst0.empty())
		{
			printf_s("layer= %d \r\n", mesh.data(mylst0.front()).layer);
			for (list<VertexHandle>::iterator lter = mylst0.begin(); lter != mylst0.end(); lter++)
			{
				find_nextlayer(*lter);
			}
			mylst0 = mylst1;
			mylst1.clear();
		}
	}
	void find_shortestway(VertexHandle v)
	{
		if (mesh.data(v).weight==0)
		{
			return;
		}
		
		float temp = 0; PMyMesh::VertexVertexIter vvter = mesh.vv_begin(v);	
		
		VertexHandle vtemp= (*vvter);
		for (temp=mesh.data(vvter).weight;vvter.is_valid();vvter++)
		{
			if (mesh.data(vvter).weight<temp)
			{
				temp = mesh.data(vvter).weight;
				vtemp = (*vvter);
			}

		}
		mesh.data(mesh.edge_handle(halfedgeh(v, vtemp, mesh))).sharp_ = true;
		find_shortestway(vtemp);
	
	}
	void setkey()
	{
		int vh_ti = 0;
		for (PMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{   //cout << (*vter).idx() << endl;
			double value = mesh.data(vter).f_;
			bool temp1 = true; int temp2 = 0;
			PMyMesh::VertexVertexIter vvter = mesh.vv_iter(*vter);
			if (value<mesh.data(vvter).f_)
			{
				temp1 = false;
			}
			else if (value>mesh.data(vvter).f_)
			{
				temp1 = true;
			}
			vvter++;
			for (; vvter.is_valid(); vvter++)
			{
				if (temp1 == true && value <= mesh.data(vvter).f_)
				{
					temp1 = false;
					temp2++;
				}
				else if (temp1 == false && value >= mesh.data(vvter).f_)
				{
					temp1 = true;
					temp2++;
				}
			}
			if (temp1 == true && value <= mesh.data(vvter).f_)
			{
				temp1 = false;
				temp2++;
			}
			else if (temp1 == false && value >= mesh.data(vvter).f_)
			{
				temp1 = true;
				temp2++;
			}

			if (temp2 == 0)
			{
				// my++;
				if (temp1 == true)
				{
					mesh.data(vter).key_ = (1);
					//strncpy(mesh.data(vter).key(),"max",3);
				}
				else
				{
					mesh.data(vter).key_ = (2);
				}
			}
			else if (temp2 >= 4)
			{
				mesh.data(vter).key_ = (temp2);
				vh[vh_ti] = (*vter);
				vh_ti++;
			}

		}
	}
	void findsaddle(bool ini, PMyMesh::VertexHandle v)
	{
		cout << "下一个马鞍点为：" << (v).idx() << " key  " << mesh.data(v).key_ << endl;
		if (ini)
		{
			if (mesh.data(v).key_ == 1)
			{
				return;
			}
		}
		else {
			if (mesh.data(v).key_ == 2)
			{
				return;
			}
		}

		double valu = mesh.data(v).f_;
		//double temp = valu; 
		OpenMesh::HalfedgeHandle he = mesh.halfedge_handle(v);
		//HalfedgeHandle he1 = mesh.prev_halfedge_handle(he);
		//he1 = mesh.opposite_halfedge_handle(he1);
		//cout << "v    " << mesh.data(mesh.from_vertex_handle(he)).f() << "   v     " << mesh.data(mesh.from_vertex_handle(he1)).f() << endl;
		int scar = mesh.valence(v);
		VertexHandle va[N], va1[N]; int i = 0;
		cout << "赋值" << endl;
		for (PMyMesh::VertexVertexIter vvter = mesh.vv_iter(v); vvter.is_valid(); vvter++)
		{
			va[i] = (*vvter);
			i++;
		}
		if (i != scar)
		{
			cout << "错误" << endl;
			return;
		}
		i = 0;
		cout << "分析山脊路径" << endl;
		for (int j = scar; j<2 * scar; j++)
		{
			if (ini)
			{
				if (mesh.data(va[j%scar]).f_ > valu&&mesh.data(va[(j - 1) % scar]).f_<mesh.data(va[j%scar]).f_ && mesh.data(va[(j + 1) % scar]).f_<mesh.data(va[j%scar]).f_)
				{
					va1[i] = va[j%scar];
					i++;
				}
			}
			else
			{
				if (mesh.data(va[j%scar]).f_ < valu&&mesh.data(va[(j - 1) % scar]).f_>mesh.data(va[j%scar]).f_ && mesh.data(va[(j + 1) % scar]).f_>mesh.data(va[j%scar]).f_)
				{
					va1[i] = va[j%scar];
					i++;
				}
			}
		}
		int k = 0;
		cout << "k为   " << k << "  i为  " << i << "寻找下一个点" << endl;
		while (k<i)
		{
			cout << "while" << endl;

			mesh.data(mesh.edge_handle(halfedgeh(v, va1[k],mesh))).sharp_ = (true);
			cout << "findsaddle" << endl;
			findsaddle(ini, va1[k]);
			if (mesh.data(v).key_ <= 2)
			{
				break;
			}
			k++;
		}
	}
	void drawedge()
	{
		Eigen::RowVector3d v1, v2;
		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			if (mesh.data(*eter).sharp_)
			{
				PMyMesh::HalfedgeHandle he = mesh.s_halfedge_handle(*eter, 0);
				PMyMesh::Point p1, p2;
				p1 = mesh.point(mesh.from_vertex_handle(he)); p2 = mesh.point(mesh.to_vertex_handle(he));
				v1(0) = 1.02*p1[0]; v2(0) = 1.02*p2[0];
				v1(1) = 1.02*p1[1]; v2(1) = 1.02*p2[1];
				v1(2) = 1.02*p1[2]; v2(2) = 1.02*p2[2];
				viewer.data.add_edges(v1, v2, Eigen::RowVector3d(0.5, 1, 0));
			}
		}
	}
	bool pre_draw(igl::viewer::Viewer &viewer)
	{
		if (viewer.core.is_animating)
		{
			viewer.data.clear();
			drawedge();
		}
		viewer.core.is_animating = false;
		return false;
	}
	void cutloop(VertexHandle v)
	{
		mesh.request_face_status();
		mesh.request_edge_status();
		mesh.request_vertex_status();
		HalfedgeHandle he = mesh.halfedge_handle(v), he1;
		EdgeHandle e = mesh.edge_handle(he);
		FaceHandle f = mesh.face_handle(he), fn = f; VertexHandle vhandle[3], v1 = v, vtest; PMyMesh::Point p;
		while (!mesh.data(e).sharp_)
		{
			he = mesh.ccw_rotated_halfedge_handle(he);
			e = mesh.edge_handle(he);
		}

		if (!mesh.data(e).sharp_)
		{
			printf_s("错误\r\n");
			return;
		}
		if (mesh.from_vertex_handle(he) != v)
		{
			printf_s("错误\r\n");
			return;
		}

		int numf = mesh.n_faces();
		Eigen::VectorXd z(numf);
		z.setZero();
		Eigen::MatrixXd c;

		int times = 0;
		//vhandle[0] = mesh.to_vertex_handle(he);
		do {
			p = mesh.point(v1);
			vhandle[0] = mesh.add_vertex(p);
			if (times == 0)
			{

				vtest = vhandle[0];
				//std::cout << "vtest id : " << vtest.idx() << std::endl;
			}
			times++;
			std::cout << times << std::endl;
			do {
				//he = he1;
				f = mesh.face_handle(he);
				vhandle[1] = mesh.to_vertex_handle(he);
				vhandle[2] = mesh.to_vertex_handle(mesh.next_halfedge_handle(he));
				he = mesh.ccw_rotated_halfedge_handle(he);
				//e = mesh.edge_handle(he);

				if (!mesh.is_valid_handle(f))
				{
					printf_s("错误\r\n"); 
					return;
				}
				mesh.delete_face(f, false);
				mesh.garbage_collection();
				fn = mesh.add_face(vhandle[0], vhandle[1], vhandle[2]);
			} while (!mesh.data(mesh.edge_handle(he)).sharp_&&mesh.is_valid_handle(mesh.face_handle(he)));

			he = mesh.halfedge_handle(fn);
			while (mesh.to_vertex_handle(he) != vhandle[0])
			{
				he = mesh.next_halfedge_handle(he);
			}
			v1 = mesh.from_vertex_handle(he);
			//std::cout << times << std::endl;
			//std::cout << "v1 id : " << v1.idx() << std::endl;
		} while (mesh.from_vertex_handle(he) != vtest);
		mesh.release_edge_status();
		mesh.release_face_status();
		mesh.release_vertex_status();
	}
	void draw_scene()
	{
		
		Eigen::VectorXf z(mesh.n_vertices());
		for (PMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			z.coeffRef((*vter).idx()) = mesh.data(*vter).key_;
			//cout << (*vter).idx() << "   " << mesh.data(*vter).key_ << "   " << mesh.data(vter).f_ << endl;
		}
		Eigen::MatrixXd c;
		igl::jet(z, true, c);
		viewer.data.set_mesh(V, F);
		viewer.data.set_colors(c);
		viewer.core.align_camera_center(V, F);
	}
	void set_frame()
	{
		viewer.callback_init = [&](igl::viewer::Viewer &viewer)
		{viewer.ngui->addVariable("kk", kk);
		viewer.ngui->addVariable<bool>("boolvariable", [&](bool var) {boolvariable = var;  },
			[&]() {return boolvariable; });
		viewer.ngui->addButton("drawedge", [&]() {viewer.data.clear(); init_sharp(); findsaddle(boolvariable, vh[kk]); drawedge(); draw_scene(); });

		viewer.screen->performLayout();
		return false;
		};
	}
	Eigen::MatrixXd A, V;
	Eigen::MatrixXi F;
private:
	list<VertexHandle> mylst0,mylst1;
	PMyMesh &mesh;
	igl::viewer::Viewer &viewer;
	int kk = 0; bool boolvariable = false;
	VertexHandle vh[N];
};


