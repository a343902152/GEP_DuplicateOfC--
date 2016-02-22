import Model.BinaryTreeNode;
import Model.Individual;

import java.util.Random;

public class Main {

    /**
     * 一个基因序列中常数的个数
     */
    public static final int constantSCount=10;

    // 测试数据的row行，col列
    public static int dataRow;
    public static int dataCol;

    // 基因重组和变异几率
    public static double pRecombination;
    public static double pMutation;

    // 系数的取值范围
    public static double lowerBound;
    public static double upperBound;


    // 物种代数、每代个体数、下代个体数
    public static int GEPNum;
    public static int groupSize;
    public static int nextGroupSize;

    // 基因数(用于多基因编程)、头部长度、尾部长度、单个基因长度，基因序列总长度(基因序列包含多个基因)
    public static int numOfGenes;
    public static int headLength;
    public static int tailLength;
    public static int perGeneLength;
    public static int geneSerialLength;

    /**
     * 适应度最高的个体
     */
    public static Individual bestFit;
    public static Individual[][] oldPopulation;
    public static Individual[] newPopulation;
    public static BinaryTreeNode[] binaryTree;
    public static double[] midResult;


    public static double[][] trainingset;
    public static char[] functionset={'+','-','*','/','@','#','$','~'};

    public static void ReadData(String inFilename,String outFilename){
        // 读取文件中的数据
        GEPNum=200;
        groupSize=500;
        nextGroupSize=100;

        numOfGenes=3;
        headLength=6;
        tailLength=headLength+1;
        perGeneLength=headLength+tailLength;
        geneSerialLength=perGeneLength*(numOfGenes+1);

        pRecombination=0.4;
        pMutation=0.15;

        lowerBound=-5;
        upperBound=5;

        dataRow=6;
        dataCol=3;

        trainingset=new double[dataRow][dataCol];
        trainingset[0][0]=1;trainingset[0][1]=1;trainingset[0][2]=5;
        trainingset[0][0]=1;trainingset[0][1]=2;trainingset[0][2]=8;
        trainingset[0][0]=0.5;trainingset[0][1]=-0.5;trainingset[0][2]=3.5;
        trainingset[0][0]=0.88 ;trainingset[0][1]=7.31 ;trainingset[0][2]=57.2105;
        trainingset[0][0]=10;trainingset[0][1]=15;trainingset[0][2]=328;
        trainingset[0][0]=-7;trainingset[0][1]=-9;trainingset[0][2]=133;

    }

    public static void InitGEP(){
        bestFit=new Individual(constantSCount);

        oldPopulation=new Individual[GEPNum+1][groupSize];
        for(int i=0;i<GEPNum+1;i++){
            for(int j=0;j<groupSize;j++){
                oldPopulation[i][j]=new Individual(constantSCount);
            }
        }

        newPopulation=new Individual[nextGroupSize];
        for(int i=0;i<nextGroupSize;i++){
            newPopulation[i]=new Individual(constantSCount);
        }

        binaryTree=new BinaryTreeNode[perGeneLength];
        for(int i=0;i<perGeneLength;i++){
            binaryTree[i]=new BinaryTreeNode();
        }
        midResult=new double[numOfGenes];

//        data=new double[nColumn];
//        infixmid=new char *[m_nNumOfGenes+1];
//        for(i=0;i<=m_nNumOfGenes;i++)
//            infixmid[i]=new char[nLengthOfGene+2*m_nHeadSize];
//        infixfinal=new char[(m_nNumOfGenes+1)*(nLengthOfGene+2*m_nHeadSize)];
//        infixlength=new int[m_nNumOfGenes+1];


        // 设置可选用的基因
        String fandt="?";
        for(int i=0;i<dataCol-1;i++){
            fandt+=(char)('A'+(i));
        }
        for(int i=0;i<functionset.length;i++){
            fandt+=functionset[i];
        }
        for(int i=0;i<numOfGenes;i++){
            fandt+=(char)('1'+i);
        }
        System.out.println(fandt);

        // 初始化种群
        // 生成初代基因序列
        for(int i=0;i<groupSize;i++){
            // 生成第123个基因的序列
            for(int j=0;j<numOfGenes;j++){
                int k;
                for(k=0;k<headLength;k++){
                    int pos=RandomInt(0,functionset.length+dataCol-1);
                    oldPopulation[0][i].geneSerial+=fandt.charAt(pos);
                }
                for(;k<perGeneLength;k++){
                    int pos=RandomInt(0,dataCol-1);
                    oldPopulation[0][i].geneSerial+=fandt.charAt(pos);
                }
            }
            // 生成合成多基因的表达式
            int k;
            for(k=0;k<headLength;k++){
                int pos=RandomInt(dataCol,dataCol+functionset.length+numOfGenes-1);
                oldPopulation[0][i].geneSerial+=fandt.charAt(pos);
            }
            for(;k<perGeneLength;k++){
                int pos=RandomInt(dataCol+functionset.length,dataCol+functionset.length+numOfGenes-1);
                oldPopulation[0][i].geneSerial+=fandt.charAt(pos);
            }
            System.out.println(oldPopulation[0][i].geneSerial);

            // 初始化系数列表
            for(int j=0;j<constantSCount;j++){
                oldPopulation[0][i].constants[j]=RandomFloat(lowerBound,upperBound);
            }
        }

        // 计算初代种群的适应度
        for(int i=0;i<groupSize;i++){
            oldPopulation[0][i].fitness= GetIndividualFitness(oldPopulation[0][i]);
            System.out.println(oldPopulation[0][i].fitness);
//            TotalF=oldpop[nCurGepNum][i].fitness+TotalF;
        }
//        for(i=0;i<=m_nNumOfGenes;i++)
//            infixlength[i]=0;
    }

    /**
     * 获取某个个体的适应度
     * 需要计算该个体的基因跟训练数据的吻合程度，满分100分
     * @param individual
     * @return
     */
    public static double GetIndividualFitness(Individual individual) {
        // fixme 容易出现正负无穷大？？？？？
        double ans=0;
        for(int i=0;i<dataRow;i++){
            int j=0;
            for(j=0;j<numOfGenes;j++){
                String serial=individual.geneSerial.substring(j*perGeneLength,(j+1)*perGeneLength);
                midResult[j]=DoGetPartFitness(i,serial,individual);
            }
            String serial=individual.geneSerial.substring(j*perGeneLength,(j+1)*perGeneLength);
            ans+=100-Math.abs((Math.abs(DoGetPartFitness(i,serial,individual)-trainingset[i][dataCol-1]))
                    /trainingset[i][dataCol-1])*100;
        }
        return ans/dataRow;
    }

    private static double DoGetPartFitness(int row, String serial, Individual individual) {
        BuildTree(serial);
        return CalculateFitness(row,binaryTree[0],individual);
    }

    /**
     * 根据基因序列，建立表达树
     * @param serial
     */
    private static void BuildTree(String serial) {
        int curCount=0;
        for(int i=0;i<serial.length();i++){
            if(serial.charAt(i)=='?'){
                binaryTree[i].element=(char)('a'+curCount++);
            }else{
                binaryTree[i].element=serial.charAt(i);
            }
            binaryTree[i].lChild=null;
            binaryTree[i].rChild=null;
        }
        int parent=0;
        int child=1;
        int leaf=0;
        int notLeaf=0;
        while((leaf!=notLeaf+1)){
            if(isOperator(binaryTree[parent].element)){
                // 加减乘除运算符
                binaryTree[parent].lChild=binaryTree[child++];
                binaryTree[parent++].rChild=binaryTree[child++];
                notLeaf++;
            }else if(isSingleOperator(binaryTree[parent].element)){
                // sin等单运算符
                binaryTree[parent].lChild=binaryTree[child++];
                binaryTree[parent++].rChild=null;
            }else if(Character.isUpperCase(binaryTree[parent].element)){
                // 大写字母，表示自变量
                binaryTree[parent].lChild=null;
                binaryTree[parent++].rChild=null;
                leaf++;
            }else if(Character.isLowerCase(binaryTree[parent].element)){
                // 小写字母，表示系数
                binaryTree[parent].lChild=null;
                binaryTree[parent++].rChild=null;
                leaf++;
            }else{
                // 其他(1,2,3)表示第i个基因
                binaryTree[parent].lChild=null;
                binaryTree[parent++].rChild=null;
                leaf++;
            }
        }

    }

    private static boolean isSingleOperator(char element) {
        switch(element)
        {
            case '@':
            case '#':
            case '$':
            case '~':
                return true;
            default:
                return false;
        }
    }

    private static boolean isOperator(char element) {
        switch(element)
        {
            case '+':
            case '-':
            case '*':
            case '/':
                return true;
            default:
                return false;
        }
    }

    /**
     * 实际根据表达式(树)得到的计算结果
     * @param row
     * @param node
     * @return
     */
    private static double CalculateFitness(int row, BinaryTreeNode node, Individual individual) {
        if(isOperator(node.element)){
            switch (node.element){
                case '+':
                    return CalculateFitness(row,node.lChild,individual)+
                                CalculateFitness(row,node.rChild,individual);
                case '-':
                    return CalculateFitness(row,node.lChild,individual)-
                                CalculateFitness(row,node.rChild,individual);
                case '*':
                    return CalculateFitness(row,node.lChild,individual)*
                                CalculateFitness(row,node.rChild,individual);
                case '/':
                    if(CalculateFitness(row,node.rChild,individual)==0){
                        return 0;
                    }else{
                        return CalculateFitness(row,node.lChild,individual)/
                                    CalculateFitness(row,node.rChild,individual);
                    }
            }
        }else if(isSingleOperator(node.element)) {
            switch (node.element) {
                case '@':
                    // sqrt
                    double ans = CalculateFitness(row, node.lChild, individual);
                    if (ans <= 0) {
                        return 0;
                    } else {
                        return Math.sqrt(ans);
                    }
                case '#':
                    // sin
                    return Math.sin(CalculateFitness(row, node.lChild, individual));
                case '$':
                    // cos
                    return Math.cos(CalculateFitness(row, node.lChild, individual));
                case '~':
                    // ln
                    double ansLog = CalculateFitness(row, node.lChild, individual);
                    if (ansLog <= 0) {
                        return 0;
                    } else {
                        return Math.log(ansLog);
                    }
            }
        }else if(Character.isUpperCase(node.element)){
            return trainingset[row][(int)(node.element-'A')];
        }else if(Character.isLowerCase(node.element)){
            return individual.constants[(int)(node.element-'a')%constantSCount];
        }else{
            return midResult[node.element-'1'];
        }
        return 0;
    }

    /**
     * 开始演化
     */
    public static void Generation(){

    }

    /**
     * 每次演化完成后，显示一下当前结果
     */
    public static void ShowResult(){

    }

    public static void main(String[] args) {
        String inFilename="D:\\JAVA_project\\GEP_DuplicateOfC++\\inputData\\3.txt";
        String outFilename="D:\\JAVA_project\\GEP_DuplicateOfC++\\inputData\\res.txt";
        ReadData(inFilename,outFilename);

        InitGEP();

        for(int i=0;i<GEPNum;i++){
            Generation();
            ShowResult();
        }
        System.out.println("finish");
    }

    public static int RandomInt(int left,int right)
    {
        Random random=new Random();
        int ans=random.nextInt((right-left)+1)+left;
        return ans;
    }

    public static double RandomFloat(double left,double right)
    {
        Random random=new Random();
        double ans=random.nextDouble()*(right-left)+left;
        return ans;
    }
}
