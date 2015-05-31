package bpcobar;

import java.math.BigInteger;
import java.io.*;
import java.util.*;
import java.util.Map.Entry;

class BPCobarMain
{
    static final int P = 2;
    static final boolean DEBUG = false;

    public static void main(String[] args) throws IOException
    {
        BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

        while(true) {
            System.out.print("bpcobar> ");
            String line = in.readLine();
            line = line.trim();
            if(line.length() == 0)
                break;
            
            RSet<BPCobar> initial = parse(line);
            if(initial == null) {
                System.err.println("Syntax error.");
                continue;
            }
            System.out.println("Parsed as "+initial);
            System.out.println();
            System.out.println(vToVBoundary(initial));

            if(DEBUG) {
                System.out.println();
                System.out.println("Cache stats:");
                System.out.println("diagonalCache: "+Monom.diagonalCache.size());
                System.out.println("rightUnitCache: "+Monom.rightUnitCache.size());
                System.out.println("mToVCache: "+Monom.mToVCache.size());
                System.out.println("vToMCache: "+Monom.vToMCache.size());
            }
        }
    }

    static RSet<BPCobar> vToVBoundary(RSet<BPCobar> initial)
    {

        RSet<BPCobar> inM = new RSet<BPCobar>();
        for(Entry<BPCobar,Q> e : initial.entrySet())
            inM.add(e.getKey().vToM(), e.getValue());
        if(DEBUG) System.out.println("In m: "+inM);

        RSet<BPCobar> bound = new RSet<BPCobar>();
        for(Entry<BPCobar,Q> e : inM.entrySet())
            bound.add(e.getKey().boundary(), e.getValue());
        if(DEBUG) System.out.println("Boundary in m: "+bound);

        RSet<BPCobar> inV = new RSet<BPCobar>();
        for(Entry<BPCobar,Q> e : bound.entrySet())
            inV.add(e.getKey().mToV(), e.getValue());
        if(DEBUG) System.out.println("Boundary in v: "+inV + "\n");

        return inV;
    }

    static RSet<BPCobar> parse(String line)
    {
        String[] tok = line.split("\\+");
        RSet<BPCobar> ret = new RSet<BPCobar>();
        for(String s : tok) {
            s = s.trim();
            int idx = s.indexOf('*');
            Q coeff = Q.ONE;
            if(idx > 0)
                coeff = Q.parse(s.substring(0,idx));

            BPCobar cobar = BPCobar.parse(s.substring(idx+1), true);
            if(cobar == null) return null;

            ret.add(cobar, coeff);
        }
        return ret;
    }
}

class BPCobar implements Comparable<BPCobar>
{
    Monom/*M or V*/ coeff;
    Monom/*T*/[] entries;
    boolean inV = false;

    BPCobar(Monom/*M*/ coeff, Monom/*T*/[] entries) {
        this.coeff = coeff;
        this.entries = entries;
    }

    BPCobar(Monom/*M or V*/ coeff, Monom/*T*/[] entries, boolean inV)
    {
        this(coeff, entries);
        this.inV = inV;
    }

    BPCobar(Monom/*M*/ a, Monom/*T*/ b, Monom/*T*/[] c) {
        coeff = a;
        entries = new Monom/*T*/[c.length + 1];
        entries[0] = b;
        for(int i = 0; i < c.length; i++)
            entries[i+1] = c[i];
    }

    RSet<BPCobar> boundary() 
    {
        RSet<BPCobar> ret = new RSet<BPCobar>();

        /* first coface: 1 | right-unit of coeff */
        RSet<MonomInMT> ru = coeff.rightUnit();
        for(Entry<MonomInMT,Q> e : ru.entrySet()) {
            MonomInMT mon = e.getKey();
            ret.add(new BPCobar(mon.m, mon.t, entries), e.getValue());
        }

        /* middle cofaces: use the diagonal and then reduce with right unit */
        for(int i = 0; i < entries.length; i++) {
            RSet<DiagonalEntry> des = entries[i].diagonal();
            for(Entry<DiagonalEntry,Q> dee : des.entrySet()) {
                DiagonalEntry de = dee.getKey();
                Monom/*T*/[] after = new Monom/*T*/[entries.length - i + 1];
                after[0] = de.a;
                after[1] = de.b;
                for(int j = i+1; j < entries.length; j++)
                    after[j-i+1] = entries[j];
                ret.add(normalizeCobar(coeff, Arrays.copyOf(entries,i), de.coeff, after), Q.sign(i+1).times(dee.getValue()));
            }
        }

        /* last coface: 1 | 1 */
        ret.add(extendByOne(), Q.sign(entries.length+1)); 

        return ret;
    }

    static RSet<BPCobar> normalizeCobar(Monom/*M*/ coeff, Monom/*T*/[] before, Monom/*M*/ mid, Monom/*T*/[] after)
    {
        /* base case */
        if(before.length == 0) {
            Monom/*M*/ c = coeff.times(mid);
            return new RSet<BPCobar>(new BPCobar(c,after));
        }

        RSet<BPCobar> ret = new RSet<BPCobar>();

        /* reduction */
        RSet<MonomInMT> ru = mid.rightUnit();
        Monom/*T*/[] newbefore = Arrays.copyOf(before, before.length - 1);
        Monom/*T*/ b = before[before.length - 1];
        for(Entry<MonomInMT,Q> e : ru.entrySet()) {
            MonomInMT mon = e.getKey();
            Monom/*T*/[] newafter = new Monom/*T*/[after.length + 1];
            for(int i = 0; i < after.length; i++)
                newafter[i+1] = after[i];
            newafter[0] = b.times(mon.t);
            RSet<BPCobar> red = normalizeCobar(coeff, newbefore, mon.m, newafter);

            ret.add(red, e.getValue());
        }

        return ret;
    }

    RSet<BPCobar> mToV()
    {
        if(inV) {
            System.err.println("Cobar already in v form!");
            System.exit(1);
        }

        RSet<BPCobar> ret = new RSet<BPCobar>();

        RSet<Monom/*V*/> vcoeff = coeff.mToV();
        for(Entry<Monom/*V*/,Q> ent : vcoeff.entrySet())
            ret.add(new BPCobar(ent.getKey(), entries, true), ent.getValue());

        return ret;
    }

    RSet<BPCobar> vToM()
    {
        if(! inV) {
            System.err.println("Cobar already in m form!");
            System.exit(1);
        }

        RSet<BPCobar> ret = new RSet<BPCobar>();

        RSet<Monom/*M*/> mcoeff = coeff.vToM();
        for(Entry<Monom/*M*/,Q> ent : mcoeff.entrySet())
            ret.add(new BPCobar(ent.getKey(), entries, false), ent.getValue());

        return ret;
    }

    BPCobar extendByOne() {
        Monom/*T*/[] newentries = Arrays.copyOf(entries, entries.length + 1);
        newentries[entries.length] = Monom.ONE;
        return new BPCobar(coeff, newentries);
    }

    @Override public String toString()
    {
        String ret = "";
        if(! coeff.isOne())
            ret += coeff.toString(inV ? 'v' : 'm') + " ";
        ret += "[ ";

        boolean first = true;
        for(Monom e : entries) {
            if(!first) ret += " | ";
            first = false;
            ret += e.toString('t');
        }
        ret += " ]";

        return ret;
    }

    @Override public int compareTo(BPCobar o)
    {
        int c;
        c = entries.length - o.entries.length;
        if(c != 0) return c;
        c = coeff.compareTo(o.coeff);
        if(c != 0) return c;
        for(int i = 0; i < entries.length; i++) {
            c = entries[i].compareTo(o.entries[i]);
            if(c != 0) return c;
        }
        return 0;
    }

    static BPCobar parse(String s, boolean isV)
    {
        s = s.trim();

        Monom coeff = Monom.ONE;
        int idx = s.indexOf('[');
        if(idx > 0) {
            coeff = Monom.parse(s.substring(0,idx), isV ? "v" : "m");
            if(coeff == null) return null;
            s = s.substring(idx+1, s.length()-1);
        } else if(idx == 0) {
            s = s.substring(1,s.length()-1);
        } else {
            coeff = Monom.parse(s, isV ? "v" : "m");
            if(coeff == null) return null;
            return new BPCobar(coeff, new Monom[] {}, isV);
        }

        String[] tok = s.split("\\|");
        ArrayList<Monom> entries = new ArrayList<Monom>();
        for(String ent : tok) {
            ent = ent.trim();
            if(ent.length() == 0) continue;
            Monom m = Monom.parse(ent, "t");
            if(m == null) return null;
            entries.add(m);
        }

        return new BPCobar(coeff, entries.toArray(new Monom[] {}), isV);
    }
}

class Monom implements Comparable<Monom>
{
    public final static Monom ONE = new Monom(new int[] {});

    public int[] exp;

    Monom(int[] exp) {
        this.exp = exp;
    }

    static Monom singleton(int k, int e) {
        if(k == 0) return ONE;
        int[] exp = new int[k];
        exp[k-1] = e;
        return new Monom(exp);
    }

    Monom times(Monom o) {
        int[] a, b;
        if(exp.length > o.exp.length) {
            a = exp;
            b = o.exp;
        } else {
            a = o.exp;
            b = exp;
        }

        int[] c = Arrays.copyOf(a, a.length);
        for(int i = 0; i < b.length; i++)
            c[i] += b[i];
        return new Monom(c);
    }

    /* removes one degree of the highest entry, and returns the result. if the result is 1, returns identically ONE */
    Monom reduce() {
        if(exp[exp.length-1] > 1) {
            int[] newexp = Arrays.copyOf(exp, exp.length);
            newexp[exp.length-1]--;
            return new Monom(newexp);
        }

        /* find the next smallest entry */
        int i;
        for(i = exp.length-2; i >= 0 && exp[i] == 0; i--);
        if(i == -1) return ONE;
        int[] newexp = Arrays.copyOf(exp, i+1);
        return new Monom(newexp);
    }
    
    /* for monom in T */
    static Map<Monom,RSet<DiagonalEntry>> diagonalCache = new TreeMap<Monom,RSet<DiagonalEntry>>();
    RSet<DiagonalEntry> diagonal()
    {
        int n = exp.length; 
        if(n == 0)
            return new RSet<DiagonalEntry>(new DiagonalEntry(Monom.ONE, Monom.ONE, Monom.ONE));

        RSet<DiagonalEntry> ret = diagonalCache.get(this);
        if(ret != null) return ret;
        ret = new RSet<DiagonalEntry>();

        Monom/*T*/ next = reduce();

        if(next == ONE) { /* we're already a singleton; apply the singleton formula */
            for(int i = 0; i < n; i++)
                for(int j = 0; j <= n-i; j++)
                    ret.add(new DiagonalEntry(singleton(i,1), singleton(j, Q.pow(Q.P, i)), singleton(n-i-j, Q.pow(Q.P, i+j))), Q.ONE);
            for(int i = 1; i < n; i++) {
                RSet<DiagonalEntry> old = singleton(n-i,1).diagonal();
                old = DiagonalEntry.pow(old, Q.P, i);
                Monom/*M*/ mi = singleton(i,1);
                for(Entry<DiagonalEntry,Q> e : old.entrySet()) {
                    DiagonalEntry de = e.getKey();
                    ret.add(new DiagonalEntry(mi.times(de.coeff), de.a, de.b), e.getValue().times(Q.MINUSONE));
                }
            }

            diagonalCache.put(this,ret);
            return ret;
        }

        /* otherwise, peel a singleton off, compute on the factors */
        Monom/*T*/ single = singleton(exp.length,1);

        RSet<DiagonalEntry> nextret = next.diagonal();
        RSet<DiagonalEntry> singleret = single.diagonal();

        /* and multiply */
        for(Entry<DiagonalEntry,Q> e1 : nextret.entrySet())
            for(Entry<DiagonalEntry,Q> e2 : singleret.entrySet())
                ret.add(e1.getKey().times(e2.getKey()), e1.getValue().times(e2.getValue()));

        diagonalCache.put(this,ret);
        return ret;
    }
    
    /* for monom in M */
    static Map<Monom,RSet<MonomInMT>> rightUnitCache = new TreeMap<Monom,RSet<MonomInMT>>();
    RSet<MonomInMT> rightUnit()
    {
        int n = exp.length; 
        if(n == 0)
            return new RSet<MonomInMT>(new MonomInMT(Monom.ONE, Monom.ONE));

        RSet<MonomInMT> ret = rightUnitCache.get(this);
        if(ret != null) return ret;
        ret = new RSet<MonomInMT>();

        Monom/*M*/ next = reduce();

        if(next == ONE) { /* we're already a singleton; apply the singleton formula */
            for(int i = 0; i <= n; i++)
                ret.add(new MonomInMT(singleton(i,1), singleton(n-i, Q.pow(Q.P, i))), Q.ONE);

            rightUnitCache.put(this, ret);
            return ret;
        }

        /* otherwise, peel a singleton off, compute on the factors */
        Monom/*M*/ single = singleton(exp.length,1);

        RSet<MonomInMT> nextret = next.rightUnit();
        RSet<MonomInMT> singleret = single.rightUnit();

        /* and multiply */
        for(Entry<MonomInMT,Q> e1 : nextret.entrySet())
            for(Entry<MonomInMT,Q> e2 : singleret.entrySet())
                ret.add(e1.getKey().times(e2.getKey()), e1.getValue().times(e2.getValue()));

        rightUnitCache.put(this, ret);
        return ret;
    }

    static Map<Monom,RSet<Monom>> mToVCache = new TreeMap<Monom,RSet<Monom>>();
    RSet<Monom/*V*/> mToV()
    {
        int n = exp.length; 
        if(n == 0)
            return new RSet<Monom/*V*/>(Monom.ONE);
        
        RSet<Monom/*V*/> ret = mToVCache.get(this);
        if(ret != null) return ret;
        ret = new RSet<Monom/*V*/>();

        Monom/*M*/ next = reduce();

        if(next == ONE) { /* we're already a singleton; apply the singleton formula */
            ret.add(singleton(n,1), Q.ONEOVERP);
            for(int i = 1; i < n; i++) {
                RSet<Monom/*V*/> sub = singleton(i,1).mToV();
                Monom/*V*/ subsing = singleton(n-i, Q.pow(Q.P, i));
                for(Entry<Monom/*V*/,Q> sube : sub.entrySet())
                    ret.add(sube.getKey().times(subsing), sube.getValue().times(Q.ONEOVERP));
            }

            mToVCache.put(this,ret);
            return ret;
        }

        /* otherwise, peel a singleton off, compute on the factors */
        Monom/*M*/ single = singleton(exp.length,1);

        RSet<Monom/*V*/> nextret = next.mToV();
        RSet<Monom/*V*/> singleret = single.mToV();

        /* and multiply */
        for(Entry<Monom/*V*/,Q> e1 : nextret.entrySet())
            for(Entry<Monom/*V*/,Q> e2 : singleret.entrySet())
                ret.add(e1.getKey().times(e2.getKey()), e1.getValue().times(e2.getValue()));

        mToVCache.put(this,ret);
        return ret;
    }

    static Map<Monom,RSet<Monom>> vToMCache = new TreeMap<Monom,RSet<Monom>>();
    RSet<Monom/*M*/> vToM()
    {
        int n = exp.length; 
        if(n == 0)
            return new RSet<Monom/*M*/>(Monom.ONE);
        
        RSet<Monom/*M*/> ret = vToMCache.get(this);
        if(ret != null) return ret;
        ret = new RSet<Monom/*M*/>();

        Monom/*V*/ next = reduce();

        if(next == ONE) { /* we're already a singleton; apply the singleton formula */
            ret.add(singleton(n,1), Q.PASQ);
            for(int i = 1; i < n; i++) {
                RSet<Monom/*M*/> sub = singleton(n-i,Q.pow(Q.P,i)).vToM();
                Monom/*M*/ subsing = singleton(i,1);
                for(Entry<Monom/*M*/,Q> sube : sub.entrySet())
                    ret.add(sube.getKey().times(subsing), sube.getValue().times(Q.MINUSONE));
            }
            vToMCache.put(this,ret);
            return ret;
        }

        /* otherwise, peel a singleton off, compute on the factors */
        Monom/*V*/ single = singleton(exp.length,1);

        RSet<Monom/*M*/> nextret = next.vToM();
        RSet<Monom/*M*/> singleret = single.vToM();

        /* and multiply */
        for(Entry<Monom/*M*/,Q> e1 : nextret.entrySet())
            for(Entry<Monom/*M*/,Q> e2 : singleret.entrySet())
                ret.add(e1.getKey().times(e2.getKey()), e1.getValue().times(e2.getValue()));

        vToMCache.put(this,ret);
        return ret;
    }
    
    public String toString(char c)
    {
        if(exp.length == 0) return "1";

        String ret = "";
        for(int i = 0; i < exp.length; i++) {
            if(exp[i] == 0) continue;
            ret += c;
            ret += (i+1);
            if(exp[i] > 1) ret += "^" + exp[i];
        }
        return ret;
    }

    @Override public String toString() {
        return toString('x');
    }

    boolean isOne() {
        return (exp.length == 0);
    }

    @Override public int compareTo(Monom o)
    {
        int c;
        c = exp.length - o.exp.length;
        if(c != 0) return c;
        for(int i = 0; i < exp.length; i++) {
            c = exp[i] - o.exp[i];
            if(c != 0) return c;
        }
        return 0;
    }

    static Monom parse(String s, String c)
    {
        s = s.trim();
        if(s.equals("1")) return Monom.ONE;
        String[] tok = s.split(c);

        Monom ret = Monom.ONE;

        for(String t : tok) {
            t = t.trim();
            if(t.length() == 0) continue;
            int idx = t.indexOf('^');
            Monom next;

            try {
                if(idx == -1)
                    next = Monom.singleton(Integer.parseInt(t), 1);
                else
                    next = Monom.singleton(Integer.parseInt(t.substring(0,idx)), Integer.parseInt(t.substring(idx+1)));
            } catch(NumberFormatException err) {
                System.err.println("Error parsing monomial in "+c+": "+s+" --- tripped on token '"+t+"'");
                return null;
            }
            ret = ret.times(next);
        }

        return ret;
    }
}

class MonomInMT implements Comparable<MonomInMT>
{
    Monom/*M*/ m;
    Monom/*T*/ t;
    MonomInMT(Monom/*M*/ m, Monom/*T*/ t) {
        this.m = m;
        this.t = t;
    }
    MonomInMT times(MonomInMT o) {
        return new MonomInMT(m.times(o.m), t.times(o.t));
    }
    @Override public int compareTo(MonomInMT o) {
        int c = m.compareTo(o.m);
        if(c != 0) return c;
        return t.compareTo(o.t);
    }

    @Override public String toString() {
        return m.toString('m') + " " + t.toString('t');
    }
}

class DiagonalEntry implements Comparable<DiagonalEntry>
{
    Monom/*M*/ coeff;
    Monom/*T*/ a;
    Monom/*T*/ b;
    DiagonalEntry(Monom/*M*/ coeff, Monom/*T*/ a, Monom/*T*/ b) {
        this.coeff = coeff;
        this.a = a;
        this.b = b;
    }
    DiagonalEntry times(DiagonalEntry o) {
        return new DiagonalEntry(coeff.times(o.coeff), a.times(o.a), b.times(o.b));
    }
    @Override public int compareTo(DiagonalEntry o)
    {
        int c = coeff.compareTo(o.coeff);
        if(c != 0) return c;
        c = a.compareTo(o.a);
        if(c != 0) return c;
        return b.compareTo(o.b);
    }
    @Override public String toString()
    {
        return coeff.toString('m') + " " + a.toString('t') + " \u2297 " + b.toString('t');
    }

    static RSet<DiagonalEntry> pow(RSet<DiagonalEntry> in, int p, int i)
    {
        if(i == 0)
            return in;

        if(p != 2) {
            System.err.println("DiagonalEntry.pow() not supported for p != 2");
            System.exit(1);
        }

        RSet<DiagonalEntry> prev = pow(in, p, i-1);
        return times(prev, prev);
    }

    static RSet<DiagonalEntry> times(RSet<DiagonalEntry> a, RSet<DiagonalEntry> b)
    {
        RSet<DiagonalEntry> ret = new RSet<DiagonalEntry>();
        for(Entry<DiagonalEntry,Q> ea : a.entrySet())
            for(Entry<DiagonalEntry,Q> eb : b.entrySet())
                ret.add( ea.getKey().times(eb.getKey()), ea.getValue().times(eb.getValue()) );
        return ret;
    }
}

class Q
{
    final static Q ONE = new Q(1L);
    final static Q MINUSONE = new Q(-1L);
    final static int P = 2; /* XXX CONFIG */
    final static Q PASQ = new Q(P);
    final static Q ONEOVERP = new Q(1,P);

    private final static long BIG = (1L << 32);

    BigInteger n, d;

    Q(long n) {
        this(n,1L);
    }
    Q(long n, long d) {
        this.n = BigInteger.valueOf(n);
        this.d = BigInteger.valueOf(d);
    }
    Q(BigInteger n, BigInteger d) {
        this.n = n;
        this.d = d;
    }

    Q reduce() {
        BigInteger gcd = n.gcd(d);
        if(gcd.equals(BigInteger.ONE)) return this;
        return new Q(n.divide(gcd), d.divide(gcd));
    }

    Q plus(Q o) {
        return new Q(d.multiply(o.n).add(o.d.multiply(n)), d.multiply(o.d));
    }

    Q times(Q o) {
        if(o == ONE) return this;
        if(this == ONE) return o;
        return new Q(n.multiply(o.n), d.multiply(o.d));
    }

    static Q sign(int i) {
        if((i & 1) == 0)
            return Q.ONE;
        return Q.MINUSONE;
    }

    boolean isZero() {
        return (n.equals(BigInteger.ZERO));
    }
    boolean isOne() {
        return (n.equals(d));
    }

    static int pow(int a, int b) {
        if(b > 30)
            System.err.println("crazy power "+b);
        if(b == 0) return 1;
        if(a == 2) return (1<<b);
        return a * pow(a,b-1);
    }

    @Override public String toString() {
        Q red = reduce();
        if(red.d.equals(BigInteger.ONE)) return red.n.toString();
        else return n + "/" + d;
    }

    static Q parse(String s)
    {
        s = s.trim();
        int idx = s.indexOf('/');

        try {
        if(idx == -1)
            return new Q(Integer.parseInt(s));
        else
            return new Q(Integer.parseInt(s.substring(0,idx)), Integer.parseInt(s.substring(idx+1)));
        } catch(NumberFormatException err) {
            System.err.println("Error parsing rational: "+s);
            return null;
        }

    }
}

class RSet<T> extends TreeMap<T,Q>
{
    RSet() { }
    RSet(T t) { add(t,Q.ONE); }

    public void add(T t, Q q)
    {
        if(q.isZero()) return;

        Q oldq = get(t);
        if(oldq == null) {
            put(t,q);
            return;
        }

        Q newq = oldq.plus(q);
        if(newq.isZero())
            remove(t);
        else
            put(t,newq);
    }

    public void add(RSet<T> r, Q q)
    {
        for(Entry<T,Q> e : r.entrySet())
            add(e.getKey(), e.getValue().times(q));
    }

    @Override public String toString()
    {
        String ret = "";

        boolean first = true;
        for(Entry<T,Q> e : entrySet()) {
            if(!first) ret += "\n + ";
            first = false;

            if(! e.getValue().isOne())
                ret += e.getValue().toString() + " ";
            ret += e.getKey().toString();
        }

        return ret;
    }
}
