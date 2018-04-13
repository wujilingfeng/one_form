#include "myholomor.h"
Eigen::MatrixXd V_uv;
igl::viewer::Viewer viewer;
PMyMesh mesh;
Mos mos(mesh, viewer);
Oneform oneform(mesh, mos);
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	if (key == '1')
	{
		mos.pre_draw1();
		oneform.set_uv(V_uv);
		viewer.data.set_uv(V_uv);
		viewer.core.show_texture = true;
	}
	return false;
}
bool predr(igl::viewer::Viewer &viewer)
{
	mos.pre_draw1();
	
	return false;
}
int main()
{

	igl::readOFF("e:/off/eight.off", mos.V, mos.F);
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, "e:/off/eight.off", opt))
	{
		for (int i = 0; i<30; i++)
		{
			printf_s("¶ÁÈ¡Ê§°Ü\r\n");
		}
		return -1;
	}
	oneform.F = mos.F;
	oneform.V = mos.V;
	oneform.init_g2();
	oneform.init_w();
	oneform.init_weight();
	//std::cout << oneform.L << std::endl;
	//mos.init_A();
	mos.set_f();
	mos.setkey();
	oneform.set_omiga();
	//mos.findsaddle(mos.vh[2]);
	//oneform.one_omiga(mos.vh[2],2);
	int sum = 0,sum1=0;
	for (PMyMesh::HalfedgeIter hter=mesh.halfedges_begin();hter!=mesh.halfedges_end();hter++)
	{
		if (mesh.data(hter).omiga[2]==1)
		{
			sum++;
		}
		if (mesh.data(hter).omiga[2]==-1)
		{
			sum1++;
		}
	}
	
	printf_s("sum: %d sum1: %d\r\n",sum,sum1);

	//for (PMyMesh::EdgeIter eter=mesh.edges_begin();eter!=mesh.edges_end();eter++)
	//{
	//	printf_s("edge id is : %d\r\n",(*eter).idx());
	//}
	oneform.set_harmonic_oneform();
   //oneform.harmonic_oneform(2);
   
   //viewer.callback_post_draw = &predr;
	printf_s("%d\r\n",mesh.n_faces());
	oneform.set_cotangentbasis();
	oneform.set_cotangentstar();
	oneform.init_ww();
	//oneform.one_omigastar(2);
	for (PMyMesh::EdgeIter eter=mesh.edges_begin();eter!=mesh.edges_end();eter++)
	{
		//printf_s("he1:%f,he2:%f\r\n",mesh.data(mesh.s_halfedge_handle(*eter,0)).omiga[2], mesh.data(mesh.s_halfedge_handle(*eter, 1)).omiga[2]);

	}
	//oneform.set_omigastar();
	
	oneform.one_tromigastar(2);
	/*Eigen::Vector3f temp1, temp2;
	temp1.setZero();
	temp2.setZero();
	temp1(0) = 1;
	temp2(1) = 1;
	printf_s("%f\r\n",temp2.cross(temp1).coeff(2));*/
	//oneform.init_ww();
	//oneform.one_cotangentstar(2);
	//oneform.one_omigastar(2);
   // oneform.set_cotangentbasis();
   /*for (PMyMesh::HalfedgeIter hter=mesh.halfedges_begin();hter!=mesh.halfedges_end();hter++)
   {
	   printf_s("halfedge id is %d\r\n",(*hter).idx());
   }*/
   /*Eigen::Vector3d vv;
   vv.setOnes();
   vv.normalize();
   for (int i=0;i<3;i++)
   {
	   printf_s("%f\r\n",vv.coeff(i));
   }*/
  // printf_s("g_2= %d\r\n",2 + mesh.n_edges() - mesh.n_vertices() - mesh.n_faces());
	/*int sum = 0;
	for (PMyMesh::EdgeIter eter=mesh.edges_begin();eter!=mesh.edges_end();eter++)
	{
		if (mesh.data(eter).sharp_)
		{
			sum++;
		}

	}
	printf_s("sum= %d\r\n",sum);*/
	//mos.integral(mesh.vertex_handle(0),2);

	//mos.drawedge();
	//mos.draw_scene();
	//oneform.set_uv(V_uv);
	//Eigen::VectorXd z(mesh.n_vertices()); Eigen::MatrixXd c;
	//z = V_uv.col(0);
	//igl::jet(z,true,c);
	//viewer.data.set_colors(c);
	//std::cout << V_uv << std::endl;
	//viewer.data.set_uv(V_uv);
   // viewer.core.show_texture = true;
	
	//viewer.core.animation_max_fps = 5.;
	//init_edge_length();
	//VertexHandle v = mesh.vertex_handle(0), v1 = mesh.vertex_handle(80);
	/*mos.init_edge_length();
	mos.init_sharp();
	mos.set_destination(v);
	mos.find_shortestway(v1);
	mos.drawedge();
	mos.draw_scene();*/
   mos.ini_list(mesh.vertex_handle(100));
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = true;
	viewer.launch();
	
	/*int n = 0;
	scanf_s("d",&n);*/
	return 0;
}