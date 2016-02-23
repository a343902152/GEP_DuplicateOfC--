import Model.BinaryTreeNode;
import Model.ComparatorIndIndividual;
import Model.Individual;
import jdk.nashorn.internal.runtime.FunctionScope;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class Main {

    /**
     * 一个基因序列中常数的个数
     */
    public static final int constantSCount=10;
    public static final int _MAX=1000;

    // 测试数据的row行，col列
    public static int dataRow;
    public static int dataCol;

    // 基因重组和变异几率
    public static double pRecombination;
    public static double pMutation;

    // 系数的取值范围
    public static double lowerBound;
    public static double upperBound;


    // 物种代数、每代个体数、变异个体数
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


    // 训练集(输入的数据)、运算符集、基因序列的取值范围
    public static double[][] trainingset;
    public static char[] functionset={'+','-','*','/','@','~'}; // 加减乘除，'@'=sqrt，'~'=ln()
    public static String fandt;

    public static void ReadData(String inFilename,String outFilename){
        // 读取文件中的数据
        GEPNum=100;
        groupSize=1000;
        nextGroupSize=(int)(groupSize*(0.6));

        numOfGenes=2;
        headLength=6;
        tailLength=headLength+1;
        perGeneLength=headLength+tailLength;
        geneSerialLength=perGeneLength*(numOfGenes+1);

        pRecombination=0.4;
        pMutation=0.25;

        lowerBound=-5;
        upperBound=5;

        dataRow=6;
        dataCol=3;

        trainingset=new double[dataRow][dataCol];
        trainingset[0][0]=1;trainingset[0][1]=1;trainingset[0][2]=5;
        trainingset[1][0]=1;trainingset[1][1]=2;trainingset[1][2]=8;
        trainingset[2][0]=0.5;trainingset[2][1]=-0.5;trainingset[2][2]=3.5;
        trainingset[3][0]=0.88 ;trainingset[3][1]=7.31 ;trainingset[3][2]=57.2105;
        trainingset[4][0]=10;trainingset[4][1]=15;trainingset[4][2]=328;
        trainingset[5][0]=-7;trainingset[5][1]=-9;trainingset[5][2]=133;

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

        // 设置可选用的基因
        fandt="?";
        for(int i=0;i<dataCol-1;i++){
            fandt+=(char)('A'+(i));
        }
        for(int i=0;i<functionset.length;i++){
            fandt+=functionset[i];
        }
        for(int i=0;i<numOfGenes;i++){
            fandt+=(char)('1'+i);
        }

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
//            System.out.println(oldPopulation[0][i].geneSerial);

            // 初始化系数列表
            for(int j=0;j<constantSCount;j++){
                oldPopulation[0][i].constants[j]=RandomFloat(lowerBound,upperBound);
            }
        }

        // 计算初代种群的适应度
        for(int i=0;i<groupSize;i++){
            oldPopulation[0][i].fitness= GetIndividualFitness(oldPopulation[0][i]);
        }
        for(int i=0;i<nextGroupSize;i++){
            CopyIndividual(newPopulation[i],oldPopulation[0][i]);
        }
        SortPopulation(0);
    }

    /**
     * 获取某个个体的适应度
     * 需要计算该个体的基因跟训练数据的吻合程度，满分100分
     * @param individual
     * @return
     */
    public static double GetIndividualFitness(Individual individual) {
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
            case '^':
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
                        return _MAX;
                    }else{
                        return CalculateFitness(row,node.lChild,individual)/
                                    CalculateFitness(row,node.rChild,individual);
                    }
                case '^':
                    return Math.pow(CalculateFitness(row,node.lChild,individual),
                            CalculateFitness(row,node.rChild,individual));
                default:
                    return _MAX;
            }
        }else if(isSingleOperator(node.element)) {
            switch (node.element) {
                case '@':
                    // sqrt
                    double ans = CalculateFitness(row, node.lChild, individual);
                    if (ans <= 0) {
                        return _MAX;
                    } else {
                        return Math.sqrt(ans);
                    }
                case '~':
                    // ln
                    double ansLog = CalculateFitness(row, node.lChild, individual);
                    if (ansLog <= 0) {
                        return _MAX;
                    } else {
                        return Math.log(ansLog);
                    }
                default:
                    return _MAX;
            }
        }else if(Character.isUpperCase(node.element)){
            return trainingset[row][(int)(node.element-'A')];
        }else if(Character.isLowerCase(node.element)){
            return individual.constants[(int)(node.element-'a')%constantSCount];
        }else{
            return midResult[node.element-'1'];
        }
    }

    /**
     * 开始演化
     * @param generationCount 表示演化到第几代
     */
    public static void Generation(int generationCount) {
        double totalFitness = 0;
        for (int i = 0; i < groupSize; i++) {
            totalFitness += oldPopulation[generationCount][i].fitness;
        }
        for (int i = 0; i < nextGroupSize / 2; i++) {
            // 选定基因交换的起点和终点
            int exchangePos1 = RandomInt(0, geneSerialLength - 1);
            int exchangePos2 = RandomInt(0, geneSerialLength - 1);
            if (exchangePos1 > exchangePos2) {
                int t = exchangePos1;
                exchangePos1 = exchangePos2;
                exchangePos2 = t;
            }
            // 选择两个不同的个体，进行基因重组
            int individual1 = SelectIndividual(generationCount, totalFitness);
            int individual2 = SelectIndividual(generationCount, totalFitness);
            while (individual1 == individual2) {
                individual1 = SelectIndividual(generationCount, totalFitness);
                individual2= SelectIndividual(generationCount,totalFitness);
            }
            for (int j = 0; j < geneSerialLength; j++) {
                if (j >= exchangePos1 && j <= exchangePos2) {
                    newPopulation[2 * i].geneSerial = ReplaceCharAt(newPopulation[2 * i].geneSerial, j,
                            oldPopulation[generationCount][individual2].geneSerial.charAt(j));
                    newPopulation[2 * i + 1].geneSerial = ReplaceCharAt(newPopulation[2 * i + 1].geneSerial, j,
                            oldPopulation[generationCount][individual1].geneSerial.charAt(j));
                } else {
                    newPopulation[2 * i].geneSerial = ReplaceCharAt(newPopulation[2 * i].geneSerial, j,
                            oldPopulation[generationCount][individual1].geneSerial.charAt(j));
                    newPopulation[2 * i + 1].geneSerial = ReplaceCharAt(newPopulation[2 * i + 1].geneSerial, j,
                            oldPopulation[generationCount][individual2].geneSerial.charAt(j));
                }
            }
            for (int j = 0; j < constantSCount; j++) {
                newPopulation[2 * i].constants[j] = oldPopulation[generationCount][individual1].constants[j];
                newPopulation[2 * i + 1].constants[j] = oldPopulation[generationCount][individual2].constants[j];
            }

            // 关于变异
            while (HadSameSerial(generationCount, newPopulation[2 * i].geneSerial)) {
                newPopulation[2 * i].geneSerial = Mutation(newPopulation[2 * i].geneSerial);
            }
            while (HadSameSerial(generationCount, newPopulation[2 * i + 1].geneSerial)) {
                newPopulation[2 * i + 1].geneSerial = Mutation(newPopulation[2 * i + 1].geneSerial);
            }
            newPopulation[2 * i].fitness = GetIndividualFitness(newPopulation[2 * i]);
            newPopulation[2 * i + 1].fitness = GetIndividualFitness(newPopulation[2 * i + 1]);
        }
    }

    private static String Mutation(String str) {
        int mutationPos;
        double pm=0;
        int h;
        for(h=0;h<perGeneLength*numOfGenes;h++){
            pm=RandomFloat(0,1);
            if(pm<pMutation){
                if(h%perGeneLength<headLength){
                    mutationPos=RandomInt(0,functionset.length+dataCol-1);
                }else{
                    mutationPos=RandomInt(0,dataCol-1);
                }
                str=ReplaceCharAt(str,h,fandt.charAt(mutationPos));
            }
        }
        for(;h<str.length();h++){
            pm=RandomFloat(0,1);
            if(pm<pMutation){
                if(h%perGeneLength<headLength){
                    mutationPos=RandomInt(dataCol,dataCol+ functionset.length+numOfGenes-1);
                }else{
                    mutationPos=RandomInt(functionset.length+dataCol,dataCol+ functionset.length+numOfGenes-1);
                }
                str=ReplaceCharAt(str,h,fandt.charAt(mutationPos));
            }
        }
        return str;
    }

    private static boolean HadSameSerial(int generationCount, String str) {
        boolean hadSame=false;
        for(int j=0;j<groupSize;j++)
        {
            if(oldPopulation[generationCount][j].geneSerial.equals(str)){
                hadSame=true;
                break;
            }
        }
        return hadSame;
    }

    private static String ReplaceCharAt(String str,int index,char c){
        StringBuilder stringBuilder=new StringBuilder(str);
        stringBuilder.deleteCharAt(index);
        stringBuilder.insert(index,c);
        return stringBuilder.toString();
    }

    /**
     * 根据轮盘赌算法，选择一个合适的个体
     * @param generationCount
     * @param total
     * @return
     */
    private static int SelectIndividual(int generationCount, double total) {
        double sum=0;
        double pick=RandomFloat(0,1);
        int i;
        if(total>0){
            for(i=0;(sum<pick)&&(i<groupSize);i++)
                sum+=oldPopulation[generationCount][i].fitness/total;
        }else
            i=RandomInt(1,groupSize);
        return (i-1);
    }

    /**
     * 按照fitness从大到小排序
     * @param individuals
     * @param size
     * @return
     */
    private static List<Individual> DoSort(Individual[] individuals,int size){
        List<Individual> list=new ArrayList<>();
        for(int i=0;i<size;i++){
            list.add(individuals[i]);
        }
        ComparatorIndIndividual comparator=new ComparatorIndIndividual();
        Collections.sort(list,comparator);
        return list;
    }

    /**
     * 按适应度从打到小排列种群内的个体
     */
    private static void SortPopulation(int generationCount) {
        List<Individual> listOld=DoSort(oldPopulation[generationCount],groupSize);
        for(int i=0;i<groupSize;i++){
            oldPopulation[generationCount][i]=listOld.get(i);
        }
        List<Individual> listNew=DoSort(newPopulation,nextGroupSize);
        for(int i=0;i<nextGroupSize;i++){
            newPopulation[i]=listNew.get(i);
        }

        int a=0;
        int b=0;
        for(int insertPos=0;insertPos<groupSize;insertPos++){
            try{
                if(b<nextGroupSize&&oldPopulation[generationCount][a].fitness<newPopulation[b].fitness){
                    CopyIndividual(oldPopulation[generationCount+1][insertPos],newPopulation[b++]);
                }else{
                    CopyIndividual(oldPopulation[generationCount+1][insertPos],oldPopulation[generationCount][a++]);
                }
            }catch (Exception e){
                e.printStackTrace();
            }
        }
    }

    private static void CopyIndividual(Individual a,Individual b){
        a.fitness=b.fitness;
        a.geneSerial=b.geneSerial;
        for(int i=0;i<constantSCount;i++){
            a.constants[i]=b.constants[i];
        }
    }

    /**
     * 每次演化完成后，显示一下当前结果
     */
    public static void ShowResult(int generationCount){
        System.out.println("bestfit: "+bestFit.fitness);
    }

    public static void main(String[] args) {

        String inFilename="D:\\JAVA_project\\GEP_DuplicateOfC++\\inputData\\3.txt";
        String outFilename="D:\\JAVA_project\\GEP_DuplicateOfC++\\inputData\\res.txt";
        ReadData(inFilename,outFilename);
        InitGEP();

        List<Individual> list=new ArrayList<>();
        for(int cas=0;cas<40;cas++){
            System.out.println("Case "+cas+":");
            for(int i=0;i<GEPNum;i++){
                Generation(i);
                SortPopulation(i);
                UpdateBestFit(i);
//                ShowResult(i);
            }
            list.add(oldPopulation[GEPNum-1][0]);
            System.out.println(bestFit.fitness);
            System.out.println(bestFit.geneSerial);
        }

//        for(int i=0;i<GEPNum;i++){
//            System.out.println("generation i="+i);
//            Generation(i);
//            SortPopulation(i);
//            UpdateBestFit(i);
//            ShowResult(i);
//        }
//        System.out.println("finished");
//        System.out.println(oldPopulation[GEPNum-1][0].fitness);
//        System.out.println(oldPopulation[GEPNum-1][0].geneSerial);
//        for(int i=0;i<constantSCount;i++){
//            System.out.println(i+":"+oldPopulation[GEPNum-1][0].constants[i]);
//        }
    }

    private static void UpdateBestFit(int generationCount) {
        if(oldPopulation[generationCount][0].fitness>bestFit.fitness){
            CopyIndividual(bestFit,oldPopulation[generationCount][0]);
        }
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
