let fs = require('fs');
let readline = require('readline')
const hash = new Map();
// {
//     id:1,
//     val:[],
//     label:0,1,2,3,
// }

async function fileToObjectArray(path) {
    const stream = fs.createReadStream(path);
    const rl = readline.createInterface({
        input: stream,
    });
    let val = [];
    for await (const line of rl) {
        val.push(line.split(',').slice(0, -1).map(v => parseFloat(v)));
    }
    return val;
}

async function run(path, k, callback) {
    // 获取数据
    let data = await fileToObjectArray(path);
    let labels = data.map(() => -1);
    // 生成obj对象
    let nodes = data.map((x, i) => ([{ id: i, val: x, label: null }]));

    // 不停合并
    while (nodes.length > k) {

        // 找出相邻距离最短的两个集合
        let min = Number.POSITIVE_INFINITY;
        let minSet = [];
        let minSetIndex = [null, null]

        for (let index = 0; index < nodes.length; index++) {
            const elementI = nodes[index];
            for (let j = index + 1; j < nodes.length; j++) {
                const elementJ = nodes[j];
                let distance = callback(nodes[index], nodes[j]);
                if (distance < min) {
                    min = distance;
                    // 解压
                    minSet = [...nodes[index],...nodes[j]];
                    minSetIndex[0] = index;
                    minSetIndex[1] = j;

                }
            }
        }

        // 合并对象
        if(minSetIndex[0]<minSetIndex[1]){
            nodes.splice(minSetIndex[0], 1)
            nodes.splice(minSetIndex[1]-1, 1)
        }else{
            nodes.splice(minSetIndex[0], 1)
            nodes.splice(minSetIndex[1], 1)
        }

        nodes.push(minSet)

    }
    // 打标签
    for (let i = 0; i < nodes.length; i++) {
        const node = nodes[i];
        for (let j = 0; j < node.length; j++) {
            const element = node[j];
            labels[element.id] = i;
        }
    }
    // 输出标签
    console.log(nodes.map((x)=>{
        return x.map((x)=>x.id)
    }))
}
function singlelink(a, b) {
    let min = Number.POSITIVE_INFINITY;

    for (let i = 0; i < a.length; i++) {
        const elementA = a[i];
        for (let j = 0; j < b.length; j++) {
            const elementB = b[j];
            let result = 0;

            //使用hash表
            let idString = `${a[i].id}_${b[j].id}`;
            if (hash.has(idString)) {
                result = hash.get(idString)
            } else {
                let sum = 0;
                try {
                    for (let k = 0; k < elementA.val.length; k++) {
                        sum += Math.pow(elementA.val[k] - elementB.val[k], 2);
                    } 
                } catch (error) {
                    console.log(elementA.val,"ss",elementB)
                }

                result = Math.sqrt(sum);
                hash.set(idString, result)
            }
            if (result < min) {
                min = result
            }
        }

    }
    return min
}

function completelink (a, b) {
    let max = Number.NEGATIVE_INFINITY;

    for (let i = 0; i < a.length; i++) {
        const elementA = a[i];
        for (let j = 0; j < b.length; j++) {
            const elementB = b[j];
            let result = 0;

            //使用hash表
            let idString = `${a[i].id}_${b[j].id}`;
            if (hash.has(idString)) {
                result = hash.get(idString)
            } else {
                let sum = 0;
                try {
                    for (let k = 0; k < elementA.val.length; k++) {
                        sum += Math.pow(elementA.val[k] - elementB.val[k], 2);
                    } 
                } catch (error) {
                    console.log(elementA.val,"ss",elementB)
                }

                result = Math.sqrt(sum);
                hash.set(idString, result)
            }
            if (result < max) {
                max = result
            }
        }

    }
    return max
}

function averageLink (a, b) {
    let line = 0;

    for (let i = 0; i < a.length; i++) {
        const elementA = a[i];
        for (let j = 0; j < b.length; j++) {
            const elementB = b[j];
            let result = 0;

            //使用hash表
            let idString = `${a[i].id}_${b[j].id}`;
            if (hash.has(idString)) {
                result = hash.get(idString)
            } else {
                let sum = 0;
                try {
                    for (let k = 0; k < elementA.val.length; k++) {
                        sum += Math.pow(elementA.val[k] - elementB.val[k], 2);
                    } 
                } catch (error) {
                    console.log(elementA.val,"ss",elementB)
                }

                result = Math.sqrt(sum);
            }
            line += result;
        }

    }
    return line/(a.length*b.length)
}
run('./iris.txt', 3, completelink);