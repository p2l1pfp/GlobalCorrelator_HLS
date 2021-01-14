library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity mu_regionizer is
    generic(
            ETA_CENTER : integer 
    );
    port(
            ap_clk : IN STD_LOGIC;
            ap_rst : IN STD_LOGIC;
            ap_start : IN STD_LOGIC;
            ap_done : OUT STD_LOGIC;
            ap_idle : OUT STD_LOGIC;
            ap_ready : OUT STD_LOGIC;
            newevent : IN STD_LOGIC;
            mu_in_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            mu_in_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            mu_out_0_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_1_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_2_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_3_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_4_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_5_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_6_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_7_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_8_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            mu_out_valid_0 : OUT STD_LOGIC;
            mu_out_valid_1 : OUT STD_LOGIC;
            mu_out_valid_2 : OUT STD_LOGIC;
            mu_out_valid_3 : OUT STD_LOGIC;
            mu_out_valid_4 : OUT STD_LOGIC;
            mu_out_valid_5 : OUT STD_LOGIC;
            mu_out_valid_6 : OUT STD_LOGIC;
            mu_out_valid_7 : OUT STD_LOGIC;
            mu_out_valid_8 : OUT STD_LOGIC;
            newevent_out : OUT STD_LOGIC

    );
end mu_regionizer;

architecture Behavioral of mu_regionizer is
    constant NREGIONS  : natural := NPFREGIONS;
    constant NALLFIFOS : natural := NPFREGIONS*NMUFIBERS;

    signal links_in :       glbparticles(NMUFIBERS-1 downto 0) := (others => null_glbparticle);
    signal fifo_in :        particles(NALLFIFOS-1 downto 0);
    signal fifo_in_write :  std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');
    signal fifo_in_roll  :  std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');

    signal fifo_out :         particles(NALLFIFOS-1 downto 0);
    signal fifo_out_valid :   std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');
    signal fifo_out_full:     std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');
    signal fifo_out_roll:     std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');

    signal merged_out :        particles(NREGIONS-1 downto 0);
    signal merged_out_valid :  std_logic_vector(NREGIONS-1 downto 0) := (others => '0');
    signal merged_out_roll:    std_logic_vector(NREGIONS-1 downto 0) := (others => '0');

begin

    router : entity work.mu_router 
                generic map(ETA_CENTER => ETA_CENTER)
                port map(ap_clk => ap_clk, 
                             enabled => ap_start,
                             newevent => newevent,
                             links_in => links_in,
                             fifo_in => fifo_in,
                             fifo_in_write => fifo_in_write,
                             fifo_in_roll  => fifo_in_roll);

    gen_fifos: for ireg in NALLFIFOS-1 downto 0 generate
        reg_buffer : entity work.rolling_fifo
                        --generic map(FIFO_INDEX => ireg+1)
                        port map(ap_clk => ap_clk, 
                                 d_in    => fifo_in(ireg),
                                 write_in  => fifo_in_write(ireg),
                                 roll   => fifo_in_roll(ireg),
                                 d_out    => fifo_out(ireg),
                                 valid_out  => fifo_out_valid(ireg),
                                 --dbg_w64 =>  fifo_dbg(ireg),
                                 full  => fifo_out_full(ireg),
                                 roll_out  => fifo_out_roll(ireg)
                             );
        end generate gen_fifos;

    gen_mergers: for imerge in NREGIONS-1 downto 0 generate
        reg_merger : entity work.fifo_merge2
                        --generic map(FIFO_INDEX => imerge+1)
                        port map(ap_clk => ap_clk, 
                                 d1_in => fifo_out(imerge*2),
                                 d2_in => fifo_out(imerge*2+1),
                                 d1_valid => fifo_out_valid(imerge*2),
                                 d2_valid => fifo_out_valid(imerge*2+1),
                                 roll     => fifo_out_roll(imerge*2),
                                 d_out      => merged_out(imerge),
                                 valid_out  => merged_out_valid(imerge),
                                 full1      => fifo_out_full(imerge*2),  
                                 full2      => fifo_out_full(imerge*2+1),
                                 --dbg_w64    => merged2_dbg(imerge),
                                 roll_out   => merged_out_roll(imerge)
                            );
        end generate gen_mergers;

    links_in( 0) <= w64_to_glbparticle(mu_in_0_V);
    links_in( 1) <= w64_to_glbparticle(mu_in_1_V);

    mu_out_0_V <= particle_to_w64(merged_out(0));
    mu_out_1_V <= particle_to_w64(merged_out(1));
    mu_out_2_V <= particle_to_w64(merged_out(2));
    mu_out_3_V <= particle_to_w64(merged_out(3));
    mu_out_4_V <= particle_to_w64(merged_out(4));
    mu_out_5_V <= particle_to_w64(merged_out(5));
    mu_out_6_V <= particle_to_w64(merged_out(6));
    mu_out_7_V <= particle_to_w64(merged_out(7));
    mu_out_8_V <= particle_to_w64(merged_out(8));
    mu_out_valid_0 <= merged_out_valid(0);
    mu_out_valid_1 <= merged_out_valid(1);
    mu_out_valid_2 <= merged_out_valid(2);
    mu_out_valid_3 <= merged_out_valid(3);
    mu_out_valid_4 <= merged_out_valid(4);
    mu_out_valid_5 <= merged_out_valid(5);
    mu_out_valid_6 <= merged_out_valid(6);
    mu_out_valid_7 <= merged_out_valid(7);
    mu_out_valid_8 <= merged_out_valid(8);

    newevent_out <= merged_out_roll(0);

end Behavioral;
